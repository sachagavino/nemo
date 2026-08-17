! ===========================================================================
! tests/coag_ice_ic_check.f90
!
! Gate (d) for the ice initial-condition machinery, checked at the INIT state
! (before any integration, so grain populations are exactly the IC ones):
!   round-trip : sum_k X_k = X_total for each base ice species (from abundances.in);
!   shape      : the area-weighted default distributes X_k proportional to n_k a_k^2,
!                so X_k/(n_k a_k^2) is CONSTANT across bins (not merely conserving the
!                total). This is the check that the *shape* is right, not just the sum.
!
! Run (no surface_ice_distribution.in present => area-weighted default) via
! tests/coag_ice_ic.sh. Expected RESULT: PASS.
! ===========================================================================
program coag_ice_ic_check
  use global_variables
  use gasgrain
  use dust_evolution
  implicit none
  integer :: b, k, idx, nbug
  character(len=11) :: base
  real(double_precision) :: tot, Xk, nk, ak, rmean, spread
  real(double_precision), allocatable :: r(:)

  call init_gasgrain()
  nbug = 0
  allocate(r(nb_grains))

  if (nb_ice_ic < 1) then
    write(*,'(a)') ' RESULT: FAIL (no deferred base-ice total -- abundances.in should give e.g. JCO=...)'
    call exit(1)
  endif

  write(*,'(a)') '=================================================================='
  write(*,'(a)') ' Ice-IC gate (d): area-weighted round-trip + shape (init state)'
  write(*,'(a)') '=================================================================='
  do b = 1, nb_ice_ic
    base = ice_ic_base_names(b)
    tot  = 0.d0
    do k = 1, nb_grains
      idx = species_index(ice_name(k, base))
      Xk  = abundances(idx)
      nk  = abundances(INDGRAIN(k))
      ak  = grain_radii(k)
      r(k) = Xk / (nk * ak*ak)
      tot  = tot + Xk
    enddo
    rmean  = sum(r) / dble(nb_grains)
    spread = (maxval(r) - minval(r)) / rmean
    write(*,'(a,a,a,es20.12,a,es9.2)') ' J', trim(base), &
      ': sum_k X_k = ', tot, '   round-trip rel-err = ', abs(tot - ice_ic_totals(b))/ice_ic_totals(b)
    write(*,'(a,es9.2,a)') '   shape  X_k/(n_k a_k^2) spread = ', spread, '  (area-weighted => ~0)'
    if (abs(tot - ice_ic_totals(b))/ice_ic_totals(b) > 1.d-10) nbug = nbug + 1
    if (spread > 1.d-10) nbug = nbug + 1
  enddo

  write(*,'(a,i0)') ' bugs = ', nbug
  if (nbug == 0) then
    write(*,'(a)') ' RESULT: PASS'
  else
    write(*,'(a)') ' RESULT: FAIL'
    call exit(1)
  endif
end program coag_ice_ic_check
