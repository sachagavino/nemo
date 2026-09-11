! ===========================================================================
! tests/coag_ice_operator_check.f90
!
! 4e diagnostic: ice-transport OPERATOR net-zero per species. With chemistry on, total
! ice is not conserved (surface reactions form/destroy it), so end-to-end ice
! conservation cannot be used. Instead we verify the COAGULATION-GENERATED ice-transport
! operator conserves ice per species BY CONSTRUCTION: for every ice-transport reaction
!   J_i X + GRAIN_j(c) -> w1 J_k1 X + w2 J_k2 X + GRAIN_j(c)
! the two products carry the SAME ice base X as the reactant and w1+w2 = 1, and the
! GRAIN partner is regenerated (catalyst, weight 1). Then the operator's net
! contribution to sum_bins Y(J_k X) is identically zero for ANY state and every step --
! a static-weight invariant, strictly stronger than a runtime sample. (jac_weight_check
! separately confirms these weights reach the RHS and Jacobian.)
!
! Run via tests/coag_ice_operator_check.sh on an ice-bearing coag fixture. RESULT: PASS.
! ===========================================================================
program coag_ice_operator_check
  use global_variables
  use gasgrain
  use dust_evolution
  implicit none
  integer :: r, nb_ice_rx, nbug
  character(len=11) :: xr, x1, x2
  real(double_precision) :: wsum, worst_wsum
  real(double_precision), parameter :: tol = 1.d-14

  call init_gasgrain()
  nb_ice_rx = 0; nbug = 0; worst_wsum = 0.d0

  do r = nb_chemistry_reactions+1, nb_reactions
    if (REACTION_TYPE(r) /= COAGULATION_TYPE) cycle
    if (REACTION_COMPOUNDS_NAMES(1,r)(1:1) /= 'J') cycle    ! ice-bearing reactant => ice-transport reaction
    nb_ice_rx = nb_ice_rx + 1

    xr = ice_base_of(REACTION_COMPOUNDS_NAMES(1,r))          ! base X of the reactant J_i X
    x1 = ice_base_of(REACTION_COMPOUNDS_NAMES(4,r))          ! base of product 1
    x2 = ice_base_of(REACTION_COMPOUNDS_NAMES(5,r))          ! base of product 2

    ! products must carry the SAME ice base as the reactant (no cross-species transport)
    if (trim(x1) /= trim(xr) .or. trim(x2) /= trim(xr)) then
      write(*,'(a,i0,5a)') '  FAIL r=',r,' base mismatch: reactant J..',trim(xr), &
        ' -> J..',trim(x1),'/J..'//trim(x2)
      nbug = nbug + 1
    endif
    ! product weights must sum to 1 (ice conserved per reaction)
    wsum = REACTION_PRODUCT_WEIGHTS(4,r) + REACTION_PRODUCT_WEIGHTS(5,r)
    worst_wsum = max(worst_wsum, abs(wsum - 1.d0))
    if (abs(wsum - 1.d0) > tol) then
      write(*,'(a,i0,a,es12.4)') '  FAIL r=',r,' w1+w2 = ', wsum
      nbug = nbug + 1
    endif
    ! partner catalyst: product 6 must equal reactant 2 with weight 1 (net-zero on GRAIN)
    if (trim(REACTION_COMPOUNDS_NAMES(6,r)) /= trim(REACTION_COMPOUNDS_NAMES(2,r)) .or. &
        abs(REACTION_PRODUCT_WEIGHTS(6,r) - 1.d0) > tol) then
      write(*,'(a,i0)') '  FAIL r=',r,' partner catalyst missing/incorrect'
      nbug = nbug + 1
    endif
  enddo

  write(*,'(a)') '=================================================================='
  write(*,'(a)') ' Ice-transport operator net-zero (structural, per species)'
  write(*,'(a)') '=================================================================='
  write(*,'(a,i0)')      ' ice-transport reactions checked = ', nb_ice_rx
  write(*,'(a,es10.2)')  ' worst |w1+w2 - 1|               = ', worst_wsum
  write(*,'(a,i0)')      ' violations                      = ', nbug
  if (nb_ice_rx > 0 .and. nbug == 0) then
    write(*,'(a)') ' RESULT: PASS'
  else
    write(*,'(a)') ' RESULT: FAIL'
    call exit(1)
  endif

contains
  !> base ice name of a per-bin ice species "J"//kk//base -> base (or '' if not ice)
  function ice_base_of(nm) result(base)
    character(len=*), intent(in) :: nm
    character(len=11) :: base
    base = ''
    if (nm(1:1) == 'J') base = trim(nm(4:))
  end function ice_base_of
end program coag_ice_operator_check
