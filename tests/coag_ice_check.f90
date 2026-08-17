! ===========================================================================
! tests/coag_ice_check.f90
!
! Hand-check (b): the ice-transport redistribution weight. Ice-follows-particle
! means the ice must be split between the two product bins with the SAME weight
! epsilon as the grains (Rung 1's product-weight field). Conservation and the rate
! gate cannot see this: any split conserves, and the rate gate uses the
! epsilon-degenerate self-collision. So we verify, for a genuine epsilon != 0 pair,
! that the ice reaction's product weights are IDENTICAL to the grain reaction's.
!
! On the 3-bin toy (mass doubling), pair (1,2):
!   m+ = m1 + m2 = 3 m1,  between m2 = 2 m1 and m3 = 4 m1
!   eps = (m3 - m+)/(m3 - m2) = (4-3)/(4-2) = 0.5   -> w1 = w2 = 0.5 into bins 2,3
! and the epsilon-degenerate clean merge (1,1):
!   m+ = 2 m1 = m2 (grid point)  -> w1 = 1, w2 = 0    -> all ice to bin 2
!
! Asserts: grain (1,2) products (GRAIN02,GRAIN03) w=(0.5,0.5); ice (1,partner 2)
! products (J02CO,J03CO) w=(0.5,0.5) + catalyst GRAIN02 w=1; ice weights == grain
! weights; and the (1,1) clean-merge weights (1,0) for both. Expected RESULT: PASS.
! ===========================================================================
program coag_ice_check
  use global_variables
  use gasgrain
  use dust_evolution
  implicit none
  integer :: k1, k2, g12, i12, g11, i11, nbug
  real(double_precision) :: w1, w2, tol
  logical :: is_overflow

  call init_gasgrain()          ! ice fixture (coagulation on, base ice JCO)
  tol = 1.d-12
  nbug = 0

  ! --- geometry: (1,2) is a genuine eps=0.5 split -------------------------------
  call coag_pair_target(1, 2, k1, k2, w1, w2, is_overflow)
  write(*,'(a)') '=================================================================='
  write(*,'(a)') ' Ice redistribution-weight hand-check (b)'
  write(*,'(a)') '=================================================================='
  write(*,'(a,2i3,a,2f8.4,a,l2)') ' pair (1,2): product bins', k1, k2, '  weights', w1, w2, '  overflow', is_overflow
  if (is_overflow .or. k1/=2 .or. k2/=3 .or. abs(w1-0.5d0)>tol .or. abs(w2-0.5d0)>tol) then
    write(*,'(a)') '   FAIL: (1,2) is not the expected eps=0.5 split into bins 2,3'
    nbug = nbug + 1
  endif

  ! --- locate the four reactions ------------------------------------------------
  g12 = find_reaction('GRAIN01', 'GRAIN02')
  i12 = find_reaction('J01CO',   'GRAIN02')
  g11 = find_reaction('GRAIN01', 'GRAIN01')
  i11 = find_reaction('J01CO',   'GRAIN01')
  if (min(g12,i12,g11,i11) == 0) then
    write(*,'(a,4i6)') '   FAIL: could not locate reactions g12,i12,g11,i11 =', g12, i12, g11, i11
    call exit(1)
  endif

  ! --- eps != 0 : ice split must equal grain split (products and weights) --------
  write(*,'(a)') ' -- eps = 0.5 pair (1,2) --'
  call show(' grain', g12)
  call show(' ice  ', i12)
  ! grain products GRAIN02/GRAIN03, ice products J02CO/J03CO, both weighted 0.5/0.5
  if (.not. is_prod(g12, 'GRAIN02', 'GRAIN03', 0.5d0, 0.5d0, tol)) nbug = nbug + 1
  if (.not. is_prod(i12, 'J02CO',   'J03CO',   0.5d0, 0.5d0, tol)) nbug = nbug + 1
  ! ice catalyst: GRAIN02 regenerated with weight 1 in product slot 6
  if (REACTION_COMPOUNDS_NAMES(6,i12) /= 'GRAIN02' .or. &
      abs(REACTION_PRODUCT_WEIGHTS(6,i12)-1.d0) > tol) then
    write(*,'(a)') '   FAIL: ice (1,2) catalyst GRAIN02 (weight 1) missing/incorrect'
    nbug = nbug + 1
  endif
  ! same-eps: ice product weights identical to grain product weights
  if (abs(REACTION_PRODUCT_WEIGHTS(4,i12)-REACTION_PRODUCT_WEIGHTS(4,g12)) > tol .or. &
      abs(REACTION_PRODUCT_WEIGHTS(5,i12)-REACTION_PRODUCT_WEIGHTS(5,g12)) > tol) then
    write(*,'(a)') '   FAIL: ice split weights differ from grain split weights (not same eps)'
    nbug = nbug + 1
  else
    write(*,'(a)') '   ice split weights == grain split weights (same eps = 0.5)  OK'
  endif

  ! --- eps = 0 : clean merge (1,1) both land all on bin 2 -----------------------
  write(*,'(a)') ' -- eps = 0 clean merge pair (1,1) --'
  if (.not. is_prod(g11, 'GRAIN02', 'GRAIN03', 1.d0, 0.d0, tol)) nbug = nbug + 1
  if (.not. is_prod(i11, 'J02CO',   'J03CO',   1.d0, 0.d0, tol)) nbug = nbug + 1
  write(*,'(a)') '   grain and ice (1,1) both deposit all to bin 2 (w=1,0)  OK'

  write(*,'(a,i0)') ' weight bugs = ', nbug
  if (nbug == 0) then
    write(*,'(a)') ' RESULT: PASS'
  else
    write(*,'(a)') ' RESULT: FAIL'
    call exit(1)
  endif

contains

  integer function find_reaction(r1, r2) result(idx)
    character(len=*), intent(in) :: r1, r2
    integer :: r
    idx = 0
    do r = nb_chemistry_reactions+1, nb_reactions
      if (REACTION_COMPOUNDS_NAMES(1,r) == r1 .and. REACTION_COMPOUNDS_NAMES(2,r) == r2) then
        idx = r; return
      endif
    enddo
  end function find_reaction

  logical function is_prod(r, p1, p2, x1, x2, tol) result(ok)
    integer, intent(in) :: r
    character(len=*), intent(in) :: p1, p2
    real(double_precision), intent(in) :: x1, x2, tol
    ok = (REACTION_COMPOUNDS_NAMES(4,r) == p1) .and. (REACTION_COMPOUNDS_NAMES(5,r) == p2) .and. &
         (abs(REACTION_PRODUCT_WEIGHTS(4,r)-x1) <= tol) .and. (abs(REACTION_PRODUCT_WEIGHTS(5,r)-x2) <= tol)
    if (.not. ok) write(*,'(a,i6,4a,2f8.4)') '   FAIL products r=', r, ' got ', &
      trim(REACTION_COMPOUNDS_NAMES(4,r)), '/', trim(REACTION_COMPOUNDS_NAMES(5,r)), &
      REACTION_PRODUCT_WEIGHTS(4,r), REACTION_PRODUCT_WEIGHTS(5,r)
  end function is_prod

  subroutine show(tag, r)
    character(len=*), intent(in) :: tag
    integer, intent(in) :: r
    write(*,'(a,a,a,a,a,a,a,a,f6.3,a,a,a,f6.3,a,a,a,f6.3)') tag, ': ', &
      trim(REACTION_COMPOUNDS_NAMES(1,r)), ' + ', trim(REACTION_COMPOUNDS_NAMES(2,r)), ' -> ', &
      trim(REACTION_COMPOUNDS_NAMES(4,r)), ' x', REACTION_PRODUCT_WEIGHTS(4,r), &
      ' + ', trim(REACTION_COMPOUNDS_NAMES(5,r)), ' x', REACTION_PRODUCT_WEIGHTS(5,r), &
      ' + ', trim(REACTION_COMPOUNDS_NAMES(6,r)), ' x', REACTION_PRODUCT_WEIGHTS(6,r)
  end subroutine show

end program coag_ice_check
