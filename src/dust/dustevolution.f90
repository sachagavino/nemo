!******************************************************************************
! MODULE: dust_evolution
!******************************************************************************
!
! DESCRIPTION:
!> @brief Rung 1: pure grain coagulation, encoded as bilinear pseudo-reactions
!!        that flow through the existing chemistry RHS/Jacobian machinery.
!!
!! ENCODING (design-thread-blessed). One pseudo-reaction per UNORDERED bin pair
!! (i,j), i<=j:
!!     bin_i + bin_j  ->  eps * bin_k  +  (1-eps) * bin_{k+1}
!! carried entirely by existing arrays -- no new RHS branch:
!!   * reactant slots 1,2 = YGRAIN(i), YGRAIN(j); slot 3 blank => no_species,
!!     so the RHS reads it as a two-body reaction (RATE = k * Y_i * Y_j * n_H).
!!   * the KERNEL K_ij rides on reaction_rates; the diagonal factor 1/2 (each
!!     unordered self-pair counted once) rides there too: reaction_rates =
!!     K_ii/2 for i==j, K_ij for i/=j. This reproduces Smoluchowski exactly:
!!       dY_i/dt = -K_ij Y_i Y_j n_H  (off-diag),  -K_ii Y_i^2 n_H  (i==j),
!!     because the RHS subtracts RATE once per reactant slot (twice when i==j).
!!   * the Podolak/Brauer REDISTRIBUTION (eps, 1-eps) rides on the per-product
!!     WEIGHT field (slots 4,5). Loss is never weighted. eps is fixed at init.
!!     Mass is conserved to machine precision because
!!       eps*m_k + (1-eps)*m_{k+1} = m_i + m_j   (linear interpolation in mass).
!!
!! TOP-BIN OVERFLOW -- Rung 1 policy (c), a LOUD guarded deferral, NOT a silent
!! omission. If m_i + m_j exceeds the top representative mass m_N, the product
!! has no bin. The long-term form (a) routes it to dust_mass_sink; Rung 1 instead
!! DOES NOT GENERATE the reaction at all -- the collision simply does not occur,
!! which is mass-exact (no leak) and, crucially, keeps the constant-kernel
!! analytic gate valid: the closed-form Smoluchowski solution has no upper
!! boundary, so any overflow would invalidate the comparison. This is only sound
!! while the top bins stay unpopulated; dust_coagulation_dropped_flux() (wired
!! with the diagnostics) reports the coagulation flux being dropped and must read
!! ~0 over the validation window, so a grid/K0 misconfiguration surfaces loudly
!! instead of hiding. The skipped pairs are recorded at init for that check.
!!
!! ICE GUARD. Coagulation moves only the refractory grain-core pseudo-species.
!! Ice (J surface / K mantle species) is not transported by coagulation until
!! Rung 2; running coagulation on an ice-bearing network would orphan the ice
!! (ice whose host grains have left the bin) -- physically incoherent. We forbid
!! it outright: coagulation on a network containing any J/K species is a fatal
!! error. Rung 1 runs on the ice-free fixture and sidesteps this entirely.
!
!******************************************************************************

module dust_evolution

use iso_fortran_env
use numerical_types
use global_variables

implicit none

! Tolerances for the geometric mass grid (mass_ratio ~ 2, well-separated bins).
real(double_precision), parameter :: COAG_MASS_TOL = 1.d-9

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Bracket the collision product m_i+m_j on the (geometric) mass grid.
!! Single source of truth for the pair -> (bins, weights) map: the count pass
!! and both fill passes call this, so they cannot disagree on which pairs exist.
!!   is_overflow = .true.  : m_i+m_j > m_N  -> no product bin (Rung 1: skip pair)
!!   otherwise product = weight w1 into bin k1, weight w2 into bin k2 (k2=k1+1),
!!   with the exact top-bin case (m_i+m_j == m_N) collapsing to a single bin.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine coag_pair_target(i, j, k1, k2, w1, w2, is_overflow)
  implicit none
  integer, intent(in)  :: i, j
  integer, intent(out) :: k1, k2
  real(double_precision), intent(out) :: w1, w2
  logical, intent(out) :: is_overflow
  real(double_precision) :: m_plus
  integer :: k, k_lo

  m_plus = mass_grid(i) + mass_grid(j)

  if (m_plus > mass_grid(nb_grains) * (1.d0 + COAG_MASS_TOL)) then
    is_overflow = .true.
    k1 = 0; k2 = 0; w1 = 0.d0; w2 = 0.d0
    return
  endif
  is_overflow = .false.

  ! largest k with m_k <= m_plus (grid is monotone increasing in mass)
  k_lo = 1
  do k = 1, nb_grains
    if (mass_grid(k) <= m_plus * (1.d0 + COAG_MASS_TOL)) k_lo = k
  enddo

  if (k_lo >= nb_grains) then
    ! lands on (or numerically at) the top bin: single product, no k+1
    k1 = nb_grains; k2 = nb_grains
    w1 = 1.d0;       w2 = 0.d0
  else
    k1 = k_lo; k2 = k_lo + 1
    ! eps = (m_{k+1} - m_plus) / (m_{k+1} - m_k)  -> weight on the lower bin
    w1 = (mass_grid(k2) - m_plus) / (mass_grid(k2) - mass_grid(k1))
    w2 = 1.d0 - w1
  endif
  return
end subroutine coag_pair_target

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Count the coagulation pseudo-reactions (non-overflow unordered pairs)
!! and record the skipped overflow pairs. Called BEFORE get_array_sizes so the
!! reaction arrays are allocated large enough; sets the module-owned counter
!! nb_coagulation_reactions that get_array_sizes adds to nb_reactions.
!! Requires mass_grid (filled in get_grain_radii, which precedes get_array_sizes).
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dust_coagulation_count()
  implicit none
  integer :: i, j, k1, k2, nover, npair
  real(double_precision) :: w1, w2
  logical :: is_overflow

  nb_coagulation_reactions = 0
  nover = 0
  do i = 1, nb_grains
    do j = i, nb_grains
      call coag_pair_target(i, j, k1, k2, w1, w2, is_overflow)
      if (is_overflow) then
        nover = nover + 1
      else
        nb_coagulation_reactions = nb_coagulation_reactions + 1
      endif
    enddo
  enddo

  coag_overflow_pairs_skipped = nover
  ! Store the skipped (i,j) so the runtime dropped-flux diagnostic can sum
  ! K_ij n_i n_j over exactly these pairs and assert it stays ~0.
  if (allocated(coag_overflow_i)) deallocate(coag_overflow_i)
  if (allocated(coag_overflow_j)) deallocate(coag_overflow_j)
  allocate(coag_overflow_i(max(nover,1)), coag_overflow_j(max(nover,1)))
  coag_overflow_i = 0; coag_overflow_j = 0
  npair = 0
  do i = 1, nb_grains
    do j = i, nb_grains
      call coag_pair_target(i, j, k1, k2, w1, w2, is_overflow)
      if (is_overflow) then
        npair = npair + 1
        coag_overflow_i(npair) = i
        coag_overflow_j(npair) = j
      endif
    enddo
  enddo

  write(*,'(a,i0,a,i0,a)') ' (coagulation) generated ', nb_coagulation_reactions, &
    ' pair reactions; ', coag_overflow_pairs_skipped, ' overflow pairs skipped (Rung 1 policy c).'
  return
end subroutine dust_coagulation_count

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Fill the coagulation reaction records at the NAME level (reactant/
!! product species names, reaction type, redistribution weights). Names are
!! resolved to species indices later by set_chemical_reactants, exactly as for
!! chemistry. MUST run after read_reactions and BEFORE index_datas (its type
!! scan needs REACTION_TYPE) and set_chemical_reactants/init_relevant_reactions
!! (the Jacobian columns depend on the resolved reactant/product IDs).
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dust_coagulation_inject_static()
  implicit none
  integer :: i, j, k1, k2, slot, ns
  real(double_precision) :: w1, w2
  logical :: is_overflow

  ! (Rung 1's coag_assert_ice_free guard is lifted here: Rung 2 adds ice-transport
  !  reactions, so an ice-bearing network is now supported rather than forbidden.)
  call coag_assert_type_unused()

  ns = 0
  do i = 1, nb_grains
    do j = i, nb_grains
      call coag_pair_target(i, j, k1, k2, w1, w2, is_overflow)
      if (is_overflow) cycle          ! Rung 1 policy (c): no reaction generated
      ns = ns + 1
      slot = nb_chemistry_reactions + ns

      REACTION_COMPOUNDS_NAMES(:, slot) = ''
      REACTION_COMPOUNDS_NAMES(1, slot) = YGRAIN(i)     ! reactant 1
      REACTION_COMPOUNDS_NAMES(2, slot) = YGRAIN(j)     ! reactant 2
      ! slot 3 (reactant 3) left blank -> no_species -> two-body in the RHS
      REACTION_COMPOUNDS_NAMES(4, slot) = YGRAIN(k1)    ! product 1 (lower bin)
      REACTION_COMPOUNDS_NAMES(5, slot) = YGRAIN(k2)    ! product 2 (upper bin)

      REACTION_TYPE(slot) = COAGULATION_TYPE
      ! Unique reaction ID (high base, clear of chemistry file IDs) so the
      ! duplicate-ID preliminary test does not treat distinct coagulation pairs
      ! as clashing entries of one reaction.
      REACTION_ID(slot) = COAG_REACTION_ID_BASE + ns

      ! Podolak/Brauer split on the product weights; loss stays unweighted.
      REACTION_PRODUCT_WEIGHTS(:, slot) = 1.d0
      REACTION_PRODUCT_WEIGHTS(4, slot) = w1
      REACTION_PRODUCT_WEIGHTS(5, slot) = w2
    enddo
  enddo

  if (ns /= nb_coagulation_reactions) then
    write(error_unit,'(a,i0,a,i0)') 'Error (coagulation): filled ', ns, &
      ' reactions but counted ', nb_coagulation_reactions
    call exit(31)
  endif
  return
end subroutine dust_coagulation_inject_static

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Set the coagulation rate coefficients (the kernel, with the diagonal
!! 1/2). MUST run after index_datas -> init_reaction_rates, which would otherwise
!! leave/overwrite these slots. Rung 1 uses the constant kernel K0; the Brownian
!! kernel (temperature-dependent, recomputed per macro-step) lands at Rung 1b.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dust_coagulation_set_rates()
  implicit none
  integer :: i, j, k1, k2, slot, ns
  real(double_precision) :: w1, w2
  logical :: is_overflow

  if (trim(coagulation_kernel) == 'constant' .and. constant_kernel_k0 <= 0.d0) then
    write(error_unit,'(a)') 'Error (coagulation): constant_kernel_k0 must be > 0 for the constant kernel.'
    call exit(31)
  endif

  ns = 0
  do i = 1, nb_grains
    do j = i, nb_grains
      call coag_pair_target(i, j, k1, k2, w1, w2, is_overflow)
      if (is_overflow) cycle
      ns = ns + 1
      slot = nb_chemistry_reactions + ns
      reaction_rates(slot) = coag_kernel(i, j)
    enddo
  enddo
  return
end subroutine dust_coagulation_set_rates

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Rate coefficient for the unordered bin pair (i,j), including the 1/2
!! factor on self-pairs (each unordered pair counted once). SINGLE SOURCE OF TRUTH
!! for the kernel: dust_coagulation_set_rates writes it onto reaction_rates and the
!! dropped-flux diagnostic evaluates it on the overflow pairs, so the guard tracks
!! whatever kernel the reactions actually use. Rung 1: constant K0. Rung 1b adds the
!! Brownian branch here (selected by coagulation_kernel); nothing else changes.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coag_kernel(i, j) result(kern)
  implicit none
  integer, intent(in) :: i, j
  real(double_precision) :: kern
  ! GRAIN coagulation rate: the bare kernel with the 1/2 self-pair factor on the
  ! diagonal (the RHS double-depletes GRAIN_i in a self-reaction, so 1/2*K_ii*2 =
  ! K_ii). Ice transport must NOT use this 1/2 (its reactants J_i X and GRAIN_j are
  ! distinct, depleted once each), so it calls coag_kernel_bare directly.
  kern = coag_kernel_bare(i, j)
  if (i == j) kern = 0.5d0 * kern             ! unordered self-pair counted once
  return
end function coag_kernel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Bare coagulation kernel K_ij for the bin pair (i,j), WITHOUT the 1/2
!! self-pair factor. Single source of truth for the kernel form. Grain reactions
!! wrap this with the diagonal 1/2 (coag_kernel); ice-transport reactions use it
!! as-is (full K_ii on the diagonal). Rung 1b adds the Brownian branch here.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coag_kernel_bare(i, j) result(kern)
  implicit none
  integer, intent(in) :: i, j
  real(double_precision) :: kern
  real(double_precision) :: ai, aj, mu

  select case (trim(coagulation_kernel))
  case ('constant')
    kern = constant_kernel_k0
  case ('brownian')
    ! Free-molecular (kinetic) Brownian coagulation kernel:
    !   K_ij = pi (a_i+a_j)^2 <v_rel>,   <v_rel> = sqrt( 8 k_B T / (pi mu_ij) )
    !        = (a_i+a_j)^2 * sqrt( 8 pi k_B T / mu_ij )        [cm^3 s^-1]
    ! <v_rel> is the MEAN (not RMS) relative thermal speed of a Maxwell-Boltzmann
    ! pair -- this fixes the prefactor convention. The gas is the momentum bath, so
    ! T is the GAS temperature (grain surface temperature is a different quantity and
    ! does not enter here). The reduced mass mu_ij = m_i m_j/(m_i+m_j) keeps the full
    ! per-bin mass dependence: the single gas T only makes the temperature factor
    ! common across pairs, it does NOT make the kernel uniform. p_stick = 1.
    ai = grain_radii(i)
    aj = grain_radii(j)
    mu = mass_grid(i) * mass_grid(j) / (mass_grid(i) + mass_grid(j))
    kern = (ai + aj)**2 * sqrt(8.d0 * PI * K_B * gas_temperature / mu)
  case default
    write(error_unit,'(3a)') 'Error (coagulation): unknown coagulation_kernel "', &
      trim(coagulation_kernel), '" (v1 supports: constant, brownian).'
    call exit(31)
  end select
  return
end function coag_kernel_bare

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Runtime overflow guard. Returns the coagulation collision rate that the
!! Rung 1 top-bin policy (c) DROPS, sum over overflow pairs of kernel(i,j) n_i n_j,
!! and the TOTAL collision rate over all pairs. The caller reports dropped/total,
!! which must stay ~0 over a valid validation window; a non-negligible value means
!! mass has reached the top bins and the boundary-free analytic solution no longer
!! holds (and, once Brownian lands, that the faster kernel is stressing the grid).
!! Uses the neutral bin population, the species the coagulation reactions act on.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dust_coagulation_flux_diag(dropped, total)
  implicit none
  real(double_precision), intent(out) :: dropped, total
  integer :: i, j, k1, k2
  real(double_precision) :: w1, w2, ni, nj, f
  logical :: is_overflow

  dropped = 0.d0
  total   = 0.d0
  do i = 1, nb_grains
    ni = abundances(INDGRAIN(i))
    do j = i, nb_grains
      nj = abundances(INDGRAIN(j))
      f  = coag_kernel(i, j) * ni * nj        ! collision-rate contribution of this pair
      total = total + f
      call coag_pair_target(i, j, k1, k2, w1, w2, is_overflow)
      if (is_overflow) dropped = dropped + f
    enddo
  enddo
  return
end subroutine dust_coagulation_flux_diag

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Guards.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine coag_assert_ice_free()
  implicit none
  integer :: s
  do s = 1, nb_species
    if (species_name(s)(1:1) == 'J' .or. species_name(s)(1:1) == 'K') then
      write(error_unit,'(a)') 'Error (coagulation): the network contains ice species (J/K) but'
      write(error_unit,'(a)') '  coagulation does not transport ice until Rung 2. Running here would'
      write(error_unit,'(a)') '  orphan the ice on grains that leave their bin. Use the ice-free'
      write(error_unit,'(a)') '  fixture for Rung 1, or wait for ice-transport coupling.'
      call exit(31)
    endif
  enddo
  return
end subroutine coag_assert_ice_free

subroutine coag_assert_type_unused()
  implicit none
  integer :: r
  do r = 1, nb_chemistry_reactions
    if (REACTION_TYPE(r) == COAGULATION_TYPE) then
      write(error_unit,'(a,i0,a)') 'Error (coagulation): reaction type ', COAGULATION_TYPE, &
        ' is already used by the chemistry network; pick another COAGULATION_TYPE.'
      call exit(31)
    endif
  enddo
  return
end subroutine coag_assert_type_unused

end module dust_evolution
