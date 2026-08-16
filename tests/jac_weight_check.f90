! ===========================================================================
! tests/jac_weight_check.f90
!
! Verifies that the per-product WEIGHT field (REACTION_PRODUCT_WEIGHTS) reaches
! every product-deposit in BOTH code paths -- the RHS get_temporal_derivatives
! AND the analytic get_jacobian -- consistently. This is the test that closes
! the gap byte-identity cannot: a weight of 1.0 is inert, so equivalence at unit
! weight proves the field is harmless but NOT that it is correctly connected.
!
! Method (differential finite difference, dependent rates frozen):
!   The weight only scales product-deposit gain terms, so the weight response
!   J(w) - J(1) isolates exactly those terms and cancels every weight-independent
!   contribution -- including the fact that get_jacobian is an APPROXIMATE Jacobian
!   (it iterates relevant_reactions(:,j) and is not the exact dF/dY even at w=1).
!   With rates frozen (freeze_dependent_rates=.true.) the true weight response is
!   FD(w) - FD(1). The two must agree to central-difference truncation on every
!   entry where get_jacobian's coverage is complete and the response is resolvable.
!
! Why the filters: get_jacobian's coverage approximation makes an ABSOLUTE FD-vs-
! analytic comparison fail at w=1 on thousands of entries that have nothing to do
! with weights. We therefore only judge entries where (a) J(1)==FD(1) to 1e-6
! (coverage provably complete here) and (b) the weight response is large enough to
! sit above the FD cancellation floor. On those, a dropped weight is unambiguous.
!
! A uniform artificial state is used so every reaction is active and every reactant
! branch (1/2/3-body, each reactant slot) is exercised at once; the identity
! J == dF/dY of the weighted product terms holds in any state.
!
! Build/run via tests/jac_weight_check.sh. Expected: RESULT: PASS.
! (Negative control: drop the weight on any one get_jacobian product deposit and
!  this reports a batch of failures with a clean ~1/w factor signature.)
! ===========================================================================
program jac_weight_check
  use global_variables
  use gasgrain
  use ode_solver
  implicit none
  integer :: n, i, j, ntest, nbug
  real(double_precision) :: t, h, colmaxJ, domcut, ybase
  real(double_precision) :: dJ, dFD, J1, Jw, FD1, FDw, reldiff, relabs1, worst
  real(double_precision) :: wb_J1, wb_FD1, wb_Jw, wb_FDw
  real(double_precision), allocatable :: Y(:), Yp(:)
  real(double_precision), allocatable :: PDJ1(:), PDJw(:), F1p(:), F1m(:), Fwp(:), Fwm(:), ian(:), jan(:)
  real(double_precision) :: w4,w5,w6,w7,w8

  ! distinct off-nominal weights, one per product slot (slots 4..8 of the record)
  w4=0.5d0; w5=0.7d0; w6=0.3d0; w7=0.2d0; w8=0.9d0
  domcut=1.d-3; ybase=1.d-4

  call init_gasgrain(); call initialize_work_arrays()
  current_time=0.d0
  call get_structure_properties(time=current_time,av=visual_extinction,density=H_number_density,gas_temperature=gas_temperature)
  call get_grain_temperature(time=current_time,gastemperature=gas_temperature,av=visual_extinction,grain_temperature=dust_temperature)
  actual_gas_temp=gas_temperature; actual_dust_temp(:)=dust_temperature(:)
  actual_av=visual_extinction; actual_gas_density=H_number_density
  n=nb_species
  allocate(Y(n),Yp(n),PDJ1(n),PDJw(n),F1p(n),F1m(n),Fwp(n),Fwm(n),ian(n),jan(n))
  ian=0.d0; jan=0.d0; Y=ybase
  ! self-shielding column densities, kept finite and consistent with the uniform state
  NH=actual_av/AV_NH_ratio*ybase; NH2=NH;NN2=NH;NCO=NH;NH2O=NH;NCH=NH;NCH3=NH
  NH2CO=NH;NCO2=NH;NN2O=NH;NCH4=NH;NOH=NH;NHCO=NH;NCN=NH;NHCN=NH;NHNC=NH;NNH=NH;NNH2=NH;NNH3=NH

  t=0.d0; freeze_dependent_rates=.false.
  call get_temporal_derivatives(n,t,Y,F1p)    ! set frozen rates at the uniform state
  freeze_dependent_rates=.true.

  ntest=0; nbug=0; worst=0.d0
  do j=1,n
    h=1.d-6*ybase
    REACTION_PRODUCT_WEIGHTS(:,:)=1.d0
    call get_jacobian(n,t,Y,j,ian,jan,PDJ1)
    Yp=Y;Yp(j)=Y(j)+h;call get_temporal_derivatives(n,t,Yp,F1p)
    Yp=Y;Yp(j)=Y(j)-h;call get_temporal_derivatives(n,t,Yp,F1m)
    REACTION_PRODUCT_WEIGHTS(4,:)=w4;REACTION_PRODUCT_WEIGHTS(5,:)=w5
    REACTION_PRODUCT_WEIGHTS(6,:)=w6;REACTION_PRODUCT_WEIGHTS(7,:)=w7;REACTION_PRODUCT_WEIGHTS(8,:)=w8
    call get_jacobian(n,t,Y,j,ian,jan,PDJw)
    Yp=Y;Yp(j)=Y(j)+h;call get_temporal_derivatives(n,t,Yp,Fwp)
    Yp=Y;Yp(j)=Y(j)-h;call get_temporal_derivatives(n,t,Yp,Fwm)
    colmaxJ=0.d0
    do i=1,n
      if(abs(PDJw(i)-PDJ1(i))>colmaxJ)colmaxJ=abs(PDJw(i)-PDJ1(i))
    end do
    if(colmaxJ<1.d-30)cycle
    do i=1,n
      dJ=PDJw(i)-PDJ1(i)
      if(abs(dJ)<domcut*colmaxJ)cycle
      J1=PDJ1(i);Jw=PDJw(i);FD1=(F1p(i)-F1m(i))/(2.d0*h);FDw=(Fwp(i)-Fwm(i))/(2.d0*h);dFD=FDw-FD1
      relabs1=abs(J1-FD1)/max(abs(J1),abs(FD1),1.d-30)
      ! judge only coverage-complete (J1==FD1) and resolvable (above FD floor) entries
      if(relabs1<1.d-6 .and. abs(dJ)>1.d-3*abs(J1) .and. abs(dJ)>1.d-13 .and. abs(J1)>1.d-12)then
        ntest=ntest+1
        reldiff=abs(dJ-dFD)/abs(dJ)
        if(reldiff>1.d-2)then
          nbug=nbug+1
          if(reldiff>worst)then;worst=reldiff;wb_J1=J1;wb_FD1=FD1;wb_Jw=Jw;wb_FDw=FDw;end if
        else
          if(reldiff>worst .and. nbug==0)worst=reldiff
        end if
      end if
    end do
  end do

  write(*,'(a)')'=================================================================='
  write(*,'(a)')' Jacobian/RHS weight-plumbing test (well-covered resolvable entries)'
  write(*,'(a)')'=================================================================='
  write(*,'(a,i7)')' entries tested = ',ntest
  write(*,'(a,i7)')' weight bugs    = ',nbug
  write(*,'(a,es12.4)')' worst reldiff  = ',worst
  if(nbug>0)then
    write(*,'(a)')' worst offending entry:'
    write(*,'(a,es13.5,a,es13.5)')'   J(1)=',wb_J1,'  FD(1)=',wb_FD1
    write(*,'(a,es13.5,a,es13.5)')'   J(w)=',wb_Jw,'  FD(w)=',wb_FDw
    write(*,'(a)')' RESULT: FAIL'
    call exit(1)
  else if(ntest>0)then
    write(*,'(a)')' RESULT: PASS'
  else
    write(*,'(a)')' RESULT: VACUOUS'
    call exit(2)
  end if
end program jac_weight_check
