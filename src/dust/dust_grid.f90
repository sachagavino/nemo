!******************************************************************************
! MODULE: dust_grid
!******************************************************************************
!
! DESCRIPTION:
!> @brief The dust size grid as a derived, first-class object.
!!
!! Two concerns are kept deliberately distinct:
!!   * the GRID       -- geometry: representative masses and radii of the bins,
!!                       derived from (a_min, a_max, mass_ratio) + grain_density.
!!                       A mass grid, geometric in mass, is the natural grid for
!!                       coagulation: a k-collision lands at m_i + m_j, and on a
!!                       geometric mass grid that target is trivially bracketed.
!!   * the DISTRIBUTION-- the initial number per bin (the IC), MRN by default,
!!                       normalised by initial_dtg_mass_ratio. This POPULATES the
!!                       grid; it does not define it.
!!
!! This module owns the derived path. The tabulated path (dust_grid_table.in,
!! the renamed legacy 0D_grain_sizes.in) is read in input_output.f90 and exists
!! so the regression against nmgc-2.0 stays bit-identical; that reference grid
!! is two representative sizes spanning ~8 decades in mass and cannot be produced
!! by any mass_ratio<=2 derivation, which is exactly why the reader is kept.
!!
!! TOP-BIN SINK (a grid property, decided now): a top-bin self-collision
!! scatters mass above m_N. On a fixed grid that mass has no bin. v1 policy:
!! it accumulates in dust_mass_sink, a diagnostic that never re-collides. A
!! non-negligible value at runtime means a_max is too low. No dynamics here yet;
!! the hook is baked into the grid so coagulation does not retrofit it.
!
!******************************************************************************

module dust_grid

use iso_fortran_env
use numerical_types
use global_variables

implicit none

!> Largest per-bin mass ratio we allow on the derived path. Coagulation
!! redistribution (Podolak/Brauer) assumes a collision product falls between two
!! adjacent bins; a mass ratio much above 2 breaks that bracketing and smears
!! mass across the grid. "<= 2" with a hair of tolerance for round-off.
real(double_precision), parameter :: MASS_RATIO_MAX = 2.d0 + 1.d-9

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Number of bins of the derived grid, from (a_min, a_max, mass_ratio).
!! Called before allocation, since it sets nb_grains. Enforces mass_ratio<=2.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dust_grid_count_bins()

  implicit none
  real(double_precision) :: m_min, m_max

  if (mass_ratio.le.1.d0) then
    write(error_unit,*) 'Error (dust_grid): mass_ratio must be > 1, got ', mass_ratio
    call exit(31)
  endif
  if (mass_ratio.gt.MASS_RATIO_MAX) then
    write(error_unit,*) 'Error (dust_grid): mass_ratio = ', mass_ratio, &
      ' exceeds the coagulation-imposed limit of 2 (mass-doubling).'
    write(error_unit,*) 'A coarser grid breaks Podolak/Brauer redistribution bracketing.'
    call exit(31)
  endif
  if (a_max.le.a_min) then
    write(error_unit,*) 'Error (dust_grid): a_max must be > a_min, got a_min=', a_min, ' a_max=', a_max
    call exit(31)
  endif

  m_min = mass_of_radius(a_min)
  m_max = mass_of_radius(a_max)

  ! m_k = m_min * mass_ratio^(k-1), k=1..nb_grains, with m_{nb_grains} <= m_max.
  ! floor(...) + 1 so the top representative mass does not exceed mass_of(a_max).
  nb_grains = int(dlog(m_max/m_min) / dlog(mass_ratio)) + 1
  if (nb_grains.lt.1) nb_grains = 1

  return
end subroutine dust_grid_count_bins

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Fill mass_grid and grain_radii for the derived grid.
!! Assumes the per-bin arrays are already allocated (nb_grains known).
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dust_grid_build_derived()

  implicit none
  integer :: k
  real(double_precision) :: m_min

  m_min = mass_of_radius(a_min)

  do k=1,nb_grains
    mass_grid(k)   = m_min * mass_ratio**(k-1)
    grain_radii(k) = radius_of_mass(mass_grid(k))
  enddo

  dust_mass_sink = 0.d0

  return
end subroutine dust_grid_build_derived

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Fill mass_grid from grain_radii for the TABULATED grid path.
!! (The table gives radii; coagulation still needs the masses.)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dust_grid_masses_from_radii()

  implicit none
  integer :: k

  do k=1,nb_grains
    mass_grid(k) = mass_of_radius(grain_radii(k))
  enddo
  dust_mass_sink = 0.d0

  return
end subroutine dust_grid_masses_from_radii

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Analytic MRN initial distribution, normalised by initial_dtg_mass_ratio.
!!
!! dn/da ~ a^p (p = dust_power_law_index, default -3.5). The number in bin k is
!! the power law integrated over the bin's radius range; the whole distribution
!! is then scaled so that the total dust mass per H equals
!! initial_dtg_mass_ratio * AMU -- the same convention as the single-grain GTODN
!! of nmgc-2.0, to which this reduces for nb_grains = 1.
!!
!! Fills GTODN_0D_temp(k) = 1 / n_k, which init_gasgrain turns into the grain
!! pseudo-species abundance abundances(INDGRAIN(k)) = 1/GTODN_0D_temp(k) = n_k.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dust_ic_mrn()

  implicit none
  integer :: k
  real(double_precision) :: p, edge_ratio_r, a_lo, a_hi
  real(double_precision) :: number_weight(nb_grains), total_mass, scale

  p = dust_power_law_index

  ! On a geometric mass grid of ratio r, bin k spans masses [m_k r^-1/2, m_k r^1/2],
  ! i.e. radii [a_k r^-1/6, a_k r^1/6]. For a tabulated grid (no fixed ratio) fall
  ! back to bin-centre bracketing by neighbours.
  edge_ratio_r = mass_ratio**(1.d0/6.d0)

  do k=1,nb_grains
    if (dust_grid_source.eq.'derived') then
      a_lo = grain_radii(k) / edge_ratio_r
      a_hi = grain_radii(k) * edge_ratio_r
    else
      if (k.eq.1) then
        a_lo = grain_radii(1)
      else
        a_lo = sqrt(grain_radii(k-1)*grain_radii(k))
      endif
      if (k.eq.nb_grains) then
        a_hi = grain_radii(nb_grains)
      else
        a_hi = sqrt(grain_radii(k)*grain_radii(k+1))
      endif
    endif

    ! integral of a^p da over [a_lo, a_hi]
    if (abs(p+1.d0).lt.1.d-12) then
      number_weight(k) = dlog(a_hi/a_lo)
    else
      number_weight(k) = (a_hi**(p+1.d0) - a_lo**(p+1.d0)) / (p+1.d0)
    endif
    if (number_weight(k).lt.0.d0) number_weight(k) = -number_weight(k)
  enddo

  ! normalise so sum_k n_k m_k = initial_dtg_mass_ratio * AMU (per H)
  total_mass = 0.d0
  do k=1,nb_grains
    total_mass = total_mass + number_weight(k) * mass_grid(k)
  enddo
  scale = initial_dtg_mass_ratio * AMU / total_mass

  do k=1,nb_grains
    GTODN_0D_temp(k) = 1.d0 / (number_weight(k) * scale)
  enddo

  return
end subroutine dust_ic_mrn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Small geometry helpers. m = (4/3) pi rho a^3.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function mass_of_radius(a) result(m)
  implicit none
  real(double_precision), intent(in) :: a
  real(double_precision) :: m
  m = (4.d0/3.d0) * PI * GRAIN_DENSITY * a*a*a
end function mass_of_radius

pure function radius_of_mass(m) result(a)
  implicit none
  real(double_precision), intent(in) :: m
  real(double_precision) :: a
  a = (3.d0*m / (4.d0*PI*GRAIN_DENSITY))**(1.d0/3.d0)
end function radius_of_mass

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Write the active grid to dust_grid_active.out at full precision, in
!! the dust_grid_table.in column format. Two purposes: (1) a diagnostic of what
!! the derivation actually produced; (2) the export side of the tabulated path
!! -- feeding this file back as dust_grid_table.in reproduces the run exactly,
!! which is how a derived grid round-trips to a restart / external hand-off.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dust_grid_export()
  implicit none
  integer :: k
  open(77, file='dust_grid_active.out', status='replace')
  write(77,'(a)') '! active dust grid (dust_grid_table.in format), written at init'
  write(77,'(a)') '! grain-radius [cm]  1/abundance  Td [K]  T_CR,peak [K]  rank'
  do k=1,nb_grains
    write(77,'(es24.16e3,2x,es24.16e3,2x,es24.16e3,2x,es24.16e3,2x,i0)') &
      grain_radii(k), GTODN_0D_temp(k), grain_temp(k), CR_PEAK_GRAIN_TEMP_all(k), k
  enddo
  close(77)
  return
end subroutine dust_grid_export

end module dust_grid
