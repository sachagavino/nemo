! ===========================================================================
! tests/ice_transport_jac_sweep.f90 -- Rung 5b DIAGNOSIS: exhaustive sweep
!
! Confirms that get_jacobian's type-50 (coagulation + ice-transport) Jacobian is
! exactly the per-reaction analytic derivative of the RHS flux, for EVERY grain
! column (all bins, both charge states) and EVERY row -- with no finite
! differencing, so no cancellation noise.
!
! Part 1 (exhaustive): for each grain column J, isolate get_jacobian to type-50
!   (zero all non-type-50 reaction_rates and the Rung-5 delta store) and compare
!   its column against an INDEPENDENT hand-summed type-50 Jacobian built by
!   looping every type-50 reaction and depositing d(flux)/dY(J) by weight/sign.
!   Run in two ice states. Any dropped reaction, wrong "other reactant", wrong
!   weight or wrong sign shows up as a mismatch.
!
! Part 2 (anchors): for a spread of (row,col) entries spanning ice-transport and
!   grain-grain across bins and charges, independently finite-difference a SUBSET
!   RHS containing only the reactions that touch that row AND depend on that col
!   (no background => no cancellation) and confirm it equals the analytic. This
!   validates the analytic flux model itself (catches e.g. a missing 2x on a
!   self-collision) at representative points.
!
! Run from a prepared run dir (coagulation = 1). Read-only: touches no src.
! ===========================================================================
program ice_transport_jac_sweep
use iso_fortran_env, only: error_unit
use global_variables
use ode_solver
use gasgrain
implicit none

integer :: N, istate, ic, ncol, r, row, k
integer, allocatable :: gcols(:)
character(len=10), allocatable :: gname(:)
real(double_precision), allocatable :: Y(:), PDJ_A(:), PDJ_B(:), IANd(:), JANd(:), rates_save(:), ydot(:)
real(double_precision) :: maxabs, maxrel, a, b, rel, tol_rel
integer :: totentries, nnz, mism, mism_total, nnz_total
character(len=24) :: statename

call init_gasgrain()
N = nb_species
allocate(Y(N),PDJ_A(N),PDJ_B(N),IANd(N),JANd(N),rates_save(nb_reactions),ydot(N))
IANd=0.d0 ; JANd=0.d0 ; tol_rel = 1.d-9

if (.not.dynamic_gtodn_active()) then
  write(error_unit,'(a)') 'FAIL: run with coagulation = 1.' ; call exit(2)
endif

! all grain columns: INDGRAIN(k)^0 and INDGRAIN_MINUS(k)^-
ncol = 2*nb_grains
allocate(gcols(ncol), gname(ncol))
do k=1,nb_grains
  gcols(k)          = INDGRAIN(k)        ; gname(k)          = trim(species_name(INDGRAIN(k)))
  gcols(nb_grains+k)= INDGRAIN_MINUS(k)  ; gname(nb_grains+k)= trim(species_name(INDGRAIN_MINUS(k)))
enddo

write(*,'(a,i0,a,i0,a,i0)') 'N=',N,'  nb_grains=',nb_grains,'  grain columns=',ncol
write(*,'(a)') '========================================================================'
write(*,'(a)') 'PART 1: get_jacobian (type-50 isolated) vs independent hand per-reaction'
write(*,'(a)') '========================================================================'

mism_total = 0 ; nnz_total = 0
do istate=1,2
  if (istate.eq.1) then ; statename='ice = 1e-10' ; else ; statename='ice = 1e-8' ; endif
  ! build state: baseline abundances, floor surface ice
  Y(1:N) = abundances(1:N)
  do row=1,N
    if (SPECIES_PHASE(row).eq.1) then
      if (istate.eq.1) then ; Y(row)=max(Y(row),1.d-10) ; else ; Y(row)=max(Y(row),1.d-8) ; endif
    endif
  enddo

  call get_temporal_derivatives(N,0.d0,Y,ydot)     ! populate chemistry rates (type-50 stay at kernel)
  rates_save(1:nb_reactions) = reaction_rates(1:nb_reactions)
  ! isolate type-50 in get_jacobian: zero non-type-50 rates + the Rung-5 delta store
  do r=1,nb_reactions
    if (REACTION_TYPE(r).ne.COAGULATION_TYPE) reaction_rates(r)=0.d0
  enddo
  if (allocated(gtodn_jac_dcoef)) gtodn_jac_dcoef(1:nb_reactions)=0.d0

  maxabs=0.d0 ; maxrel=0.d0 ; mism=0 ; nnz=0
  do ic=1,N
    call get_jacobian(N,0.d0,Y,ic,IANd,JANd,PDJ_A)
    call hand_t50_column(ic, Y, PDJ_B)
    do row=1,N
      a=PDJ_A(row) ; b=PDJ_B(row)
      if (abs(a).lt.1.d-30 .and. abs(b).lt.1.d-30) cycle
      nnz=nnz+1
      rel=abs(a-b)/max(abs(b),1.d-30)
      if (abs(a-b).gt.maxabs) maxabs=abs(a-b)
      if (rel.gt.maxrel) maxrel=rel
      if (rel.gt.tol_rel .and. abs(a-b).gt.1.d-30) then
        mism=mism+1
        if (mism.le.20) write(*,'(a,a12,a,a12,2es16.6,a,es10.2)') &
          '  MISMATCH col=',species_name(ic),' row=',species_name(row), a, b, '  rel=', rel
      endif
    enddo
  enddo
  write(*,'(a,a,a,i0,a,i0,a,es9.2,a,es9.2)') ' [',trim(statename),'] nonzero type-50 entries (all columns)=',nnz, &
    '  mismatches=',mism,'  max|A-B|=',maxabs,'  max rel=',maxrel
  mism_total=mism_total+mism ; nnz_total=nnz_total+nnz
  reaction_rates(1:nb_reactions)=rates_save(1:nb_reactions)
enddo
write(*,'(a,i0,a,i0,a)') ' PART 1 TOTAL: ',nnz_total,' nonzero entries compared, ',mism_total,' mismatches'

write(*,'(/,a)') '========================================================================'
write(*,'(a)') 'PART 2: independent subset-FD anchors (analytic model vs finite diff)'
write(*,'(a)') '========================================================================'
! rebuild state C, full rates
Y(1:N)=abundances(1:N)
do row=1,N ; if (SPECIES_PHASE(row).eq.1) Y(row)=max(Y(row),1.d-10) ; enddo
call get_temporal_derivatives(N,0.d0,Y,ydot)
rates_save(1:nb_reactions)=reaction_rates(1:nb_reactions)

! anchors spanning ice-transport (bins 2/3/4) and grain-grain (bins 2/3/4), both charge cols
call anchor('J02CO   ', INDGRAIN(1))
call anchor('J03CO   ', INDGRAIN(1))
call anchor('J04CO   ', INDGRAIN(1))
call anchor('J03H2O  ', INDGRAIN(2))
call anchor('J03CH2OH', INDGRAIN_MINUS(1))
call anchor('J04CH3OH', INDGRAIN(2))
call anchor('GRAIN02 ', INDGRAIN(1))
call anchor('GRAIN03 ', INDGRAIN(1))
call anchor('GRAIN03 ', INDGRAIN(2))
call anchor('GRAIN04 ', INDGRAIN_MINUS(2))

call exit(0)

contains

  !> Independent hand-summed type-50 Jacobian column d(RHS)/dY(J), col = J.
  subroutine hand_t50_column(J, Yin, col_out)
  integer, intent(in) :: J
  real(double_precision), intent(in) :: Yin(:)
  real(double_precision), intent(out) :: col_out(:)
  integer :: rr, s, c1, c2, c3, cc
  real(double_precision) :: df
  col_out(1:N)=0.d0
  do rr=1,nb_reactions
    if (REACTION_TYPE(rr).ne.COAGULATION_TYPE) cycle
    c1=REACTION_COMPOUNDS_ID(1,rr); c2=REACTION_COMPOUNDS_ID(2,rr); c3=REACTION_COMPOUNDS_ID(3,rr)
    if (c1.lt.1.or.c1.gt.N) cycle
    df=0.d0
    ! two-body flux = rate*Y(c1)*Y(c2)*density ; d/dY(J) picks the "other" reactant
    if (c2.ge.1.and.c2.le.N) then
      if (c1.eq.J) df = df + reaction_rates(rr)*Yin(c2)*actual_gas_density
      if (c2.eq.J) df = df + reaction_rates(rr)*Yin(c1)*actual_gas_density
    else
      if (c1.eq.J) df = df + reaction_rates(rr)                 ! one-body (not expected)
    endif
    if (df.eq.0.d0) cycle
    do s=4,8
      cc=REACTION_COMPOUNDS_ID(s,rr)
      if (cc.ge.1.and.cc.le.N) col_out(cc)=col_out(cc)+REACTION_PRODUCT_WEIGHTS(s,rr)*df
    enddo
    do s=1,3
      cc=REACTION_COMPOUNDS_ID(s,rr)
      if (cc.ge.1.and.cc.le.N) col_out(cc)=col_out(cc)-df
    enddo
  enddo
  end subroutine hand_t50_column

  !> Subset-FD anchor for entry (row=name, col=J): finite-difference a RHS that
  !! contains only the type-50 reactions that touch `row` AND depend on Y(J)
  !! (no background => no cancellation), and compare to the analytic entry.
  subroutine anchor(name, J)
  character(len=*), intent(in) :: name
  integer, intent(in) :: J
  integer :: trow, rr, s, c1, c2, c3, cc
  logical :: touch, dep
  real(double_precision) :: dY, fp, fm, fd, ana
  real(double_precision), allocatable :: Yp(:), Ym(:), hp(:), hm(:), col_out(:)
  allocate(Yp(N),Ym(N),hp(N),hm(N),col_out(N))
  trow=0
  do rr=1,N ; if (trim(species_name(rr)).eq.trim(name)) trow=rr ; enddo
  if (trow.eq.0) then ; write(*,'(a,a)') '  (anchor row not found: ',trim(name),')' ; deallocate(Yp,Ym,hp,hm,col_out) ; return ; endif

  ! analytic value for this entry
  reaction_rates(1:nb_reactions)=rates_save(1:nb_reactions)
  call hand_t50_column(J, Y, col_out) ; ana = col_out(trow)

  dY = 1.d-6*max(abs(Y(J)),1.d-30)
  Yp(1:N)=Y(1:N); Yp(J)=Y(J)+dY ; Ym(1:N)=Y(1:N); Ym(J)=Y(J)-dY
  hp(1:N)=0.d0 ; hm(1:N)=0.d0
  do rr=1,nb_reactions
    if (REACTION_TYPE(rr).ne.COAGULATION_TYPE) cycle
    c1=REACTION_COMPOUNDS_ID(1,rr); c2=REACTION_COMPOUNDS_ID(2,rr); c3=REACTION_COMPOUNDS_ID(3,rr)
    if (c1.lt.1.or.c1.gt.N) cycle
    ! subset: reaction must touch trow AND depend on Y(J)
    touch=.false.
    if (c1.eq.trow.or.c2.eq.trow.or.c3.eq.trow) touch=.true.
    do s=4,8 ; if (REACTION_COMPOUNDS_ID(s,rr).eq.trow) touch=.true. ; enddo
    dep = (c1.eq.J .or. c2.eq.J .or. c3.eq.J)
    if (.not.(touch.and.dep)) cycle
    if (c2.ge.1.and.c2.le.N) then
      fp=reaction_rates(rr)*Yp(c1)*Yp(c2)*actual_gas_density
      fm=reaction_rates(rr)*Ym(c1)*Ym(c2)*actual_gas_density
    else
      fp=reaction_rates(rr)*Yp(c1) ; fm=reaction_rates(rr)*Ym(c1)
    endif
    do s=4,8
      cc=REACTION_COMPOUNDS_ID(s,rr)
      if (cc.ge.1.and.cc.le.N) then ; hp(cc)=hp(cc)+REACTION_PRODUCT_WEIGHTS(s,rr)*fp ; hm(cc)=hm(cc)+REACTION_PRODUCT_WEIGHTS(s,rr)*fm ; endif
    enddo
    do s=1,3
      cc=REACTION_COMPOUNDS_ID(s,rr)
      if (cc.ge.1.and.cc.le.N) then ; hp(cc)=hp(cc)-fp ; hm(cc)=hm(cc)-fm ; endif
    enddo
  enddo
  fd = (hp(trow)-hm(trow))/(2.d0*dY)
  write(*,'(a,a10,a,a10,a,es15.6,a,es15.6,a,es9.2)') '  entry row=',name,' col=',trim(species_name(J)), &
    '  analytic=',ana,'  subsetFD=',fd,'  rel=',abs(ana-fd)/max(abs(fd),abs(ana),1.d-30)
  deallocate(Yp,Ym,hp,hm,col_out)
  end subroutine anchor

end program ice_transport_jac_sweep
