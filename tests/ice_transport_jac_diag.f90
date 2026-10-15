! ===========================================================================
! tests/ice_transport_jac_diag.f90 -- Rung 5b DIAGNOSIS (read-only, no fix)
!
! Characterize the existing analytic Jacobian for a SINGLE type-50 ice-transport
! reaction, term by term, against the truth. For one reaction r0 whose collision
! partner is GRAIN01^0 (so column = GRAIN01^0 is a reactant of r0):
!
!   ana(row) = hand analytic d(flux_r0)/dY(GRAIN01^0), deposited +w_p to products,
!              -1 to reactants   [the value get_jacobian SHOULD produce for r0]
!   PDJ(row) = get_jacobian's actual column, with every OTHER reaction's rate
!              coefficient and the Rung-5 delta store zeroed, so PDJ is exactly
!              r0's standard-loop contribution
!   fd (row) = central difference of a hand-built r0-only RHS w.r.t. Y(GRAIN01^0)
!              [the ground truth for r0]
!
! ana vs fd validates the hand analytic; PDJ vs fd is the verdict on the existing
! loop. Rows where PDJ disagrees with fd are exactly what the loop misses.
! Run from a prepared run dir (coagulation = 1). Read-only: touches no src.
! ===========================================================================
program ice_transport_jac_diag
use iso_fortran_env, only: error_unit
use global_variables
use ode_solver
use gasgrain
implicit none

integer :: N, r, r0, col, row, p, q, cid, nlist
real(double_precision), allocatable :: Y(:), Yp(:), Ym(:), ydot0(:), ydotp(:), ydotm(:)
real(double_precision), allocatable :: PDJ(:), fd(:), ana(:), IANd(:), JANd(:), rates_save(:)
real(double_precision) :: dY, eps, rate_r0, dflux, relerr
integer :: c(8)
real(double_precision) :: wt(8)
character(len=12) :: nm

call init_gasgrain()
N = nb_species
allocate(Y(N),Yp(N),Ym(N),ydot0(N),ydotp(N),ydotm(N),PDJ(N),fd(N),ana(N),IANd(N),JANd(N))
allocate(rates_save(nb_reactions))
IANd = 0.d0 ; JANd = 0.d0 ; eps = 1.d-6

if (.not.dynamic_gtodn_active()) then
  write(error_unit,'(a)') 'FAIL: run with coagulation = 1.' ; call exit(2)
endif

! ---- (1) enumerate type-50 ice-transport reactions (reactant1 = a surface J ice)
write(*,'(a)') '=== type-50 ice-transport reactions (first 12) ==='
nlist = 0
do r=1,nb_reactions
  if (REACTION_TYPE(r).ne.COAGULATION_TYPE) cycle
  cid = REACTION_COMPOUNDS_ID(1,r)
  if (cid.lt.1 .or. cid.gt.N) cycle
  if (SPECIES_PHASE(cid).ne.1) cycle            ! ice-transport (reactant1 is surface ice); skip grain-grain coag
  nlist = nlist + 1
  if (nlist.le.12) call print_reaction(r)
enddo
write(*,'(a,i0)') 'total ice-transport reactions: ', nlist

! ---- (2) pick r0: first ice-transport whose partner (reactant2) is GRAIN01^0
r0 = 0
do r=1,nb_reactions
  if (REACTION_TYPE(r).ne.COAGULATION_TYPE) cycle
  cid = REACTION_COMPOUNDS_ID(1,r)
  if (cid.lt.1 .or. cid.gt.N) cycle
  if (SPECIES_PHASE(cid).ne.1) cycle
  if (REACTION_COMPOUNDS_ID(2,r).eq.INDGRAIN(1)) then ; r0 = r ; exit ; endif
enddo
if (r0.eq.0) then
  write(error_unit,'(a)') 'no ice-transport reaction with GRAIN01^0 as partner found.' ; call exit(4)
endif
col = INDGRAIN(1)
write(*,'(/,a)') '=== chosen reaction r0 (column = GRAIN01^0) ==='
call print_reaction(r0)

c(1:8) = REACTION_COMPOUNDS_ID(1:8,r0)
wt(1:8) = 1.d0
wt(4:8) = REACTION_PRODUCT_WEIGHTS(4:8,r0)

! ---- (3) controlled state: baseline abundances, zero all surface ice, then set
!          r0's transported ice (reactant1) nonzero so its flux is finite.
Y(1:N) = abundances(1:N)
do row=1,N
  if (SPECIES_PHASE(row).eq.1) Y(row) = 0.d0
enddo
Y(c(1)) = 1.d-8                                  ! the ice being transported, J_iX

! ---- (4) analytic column from get_jacobian, isolated to r0 alone
call get_temporal_derivatives(N, 0.d0, Y, ydot0)  ! populate reaction_rates
rate_r0 = reaction_rates(r0)
rates_save(1:nb_reactions) = reaction_rates(1:nb_reactions)
reaction_rates(1:nb_reactions) = 0.d0
reaction_rates(r0) = rate_r0
if (allocated(gtodn_jac_dcoef)) gtodn_jac_dcoef(1:nb_reactions) = 0.d0  ! silence the Rung-5 block
call get_jacobian(N, 0.d0, Y, col, IANd, JANd, PDJ)

! ---- (5) hand analytic: d(flux_r0)/dY(col). r0 is two-body (reactant1,reactant2=col),
!          so d(flux)/dY(col) = rate * Y(reactant1) * density.
dflux = rate_r0 * Y(c(1)) * actual_gas_density
ana(1:N) = 0.d0
do p=4,8
  if (c(p).ge.1 .and. c(p).le.N) ana(c(p)) = ana(c(p)) + wt(p)*dflux
enddo
do q=1,3
  if (c(q).ge.1 .and. c(q).le.N) ana(c(q)) = ana(c(q)) - dflux
enddo

! ---- (6) hand FD of the r0-only RHS w.r.t. Y(col)
dY = eps*max(abs(Y(col)),1.d-30)
Yp(1:N)=Y(1:N); Yp(col)=Y(col)+dY
Ym(1:N)=Y(1:N); Ym(col)=Y(col)-dY
call rhs_single(r0, rate_r0, c, wt, Yp, ydotp)
call rhs_single(r0, rate_r0, c, wt, Ym, ydotm)
fd(1:N) = (ydotp(1:N)-ydotm(1:N))/(2.d0*dY)

! ---- (7) report term by term over the union of r0's compounds
write(*,'(/,a)') '=== term-by-term Jacobian column (col = GRAIN01^0) for r0 alone ==='
write(*,'(a)') 'role        species        hand-ana        get_jacobian         FD(truth)     PDJ-vs-FD'
do q=1,8
  cid = c(q)
  if (cid.lt.1 .or. cid.gt.N) cycle
  if (q.le.3) then ; nm='reactant' ; else ; nm='product' ; endif
  relerr = abs(PDJ(cid)-fd(cid))/max(abs(fd(cid)),1.d-30)
  write(*,'(a10,a13,3es16.6,a)') trim(nm), species_name(cid), ana(cid), PDJ(cid), fd(cid), &
       merge('   *** MISMATCH', '        match  ', relerr.gt.1.d-3)
enddo

! ---- (8) also scan the WHOLE column for any row where get_jacobian disagrees with
!          the r0 truth (catches deposits to rows outside r0's own compound list).
write(*,'(/,a)') '=== full-column scan: rows where PDJ (r0-isolated) disagrees with r0 FD ==='
do row=1,N
  if (abs(PDJ(row)).lt.1.d-30 .and. abs(fd(row)).lt.1.d-30) cycle
  relerr = abs(PDJ(row)-fd(row))/max(abs(fd(row)),1.d-30)
  if (relerr.gt.1.d-3) write(*,'(4x,a12,a,3es16.6)') species_name(row), '  ana/PDJ/FD=', ana(row), PDJ(row), fd(row)
enddo

! ======================================================================
! (9) STATE C decomposition: all surface ice = 1e-10 (the full-run failing
!     state). Compare full get_jacobian vs full FD vs a HAND-SUMMED type-50
!     Jacobian, on column GRAIN01^0, to see if the standard loop's accumulation
!     over many ice-transport reactions equals the per-reaction truth.
! ======================================================================
block
  real(double_precision), allocatable :: ana_t50(:)
  real(double_precision) :: df
  integer :: rr, s
  allocate(ana_t50(N))

  ! Restore the full rate vector: step (4) zeroed all but r0, and the RHS only
  ! recomputes chemistry types (not the constant-kernel type-50), so without this
  ! the type-50 rates would stay zeroed. Chemistry rates are recomputed by the
  ! get_temporal_derivatives calls below; type-50 stays at its kernel value here.
  reaction_rates(1:nb_reactions) = rates_save(1:nb_reactions)

  Y(1:N) = abundances(1:N)
  do row=1,N
    if (SPECIES_PHASE(row).eq.1 .and. Y(row).lt.1.d-10) Y(row)=1.d-10
  enddo

  ! full analytic column (real reaction_rates; Rung-5 block active)
  call get_temporal_derivatives(N,0.d0,Y,ydot0)
  call get_jacobian(N,0.d0,Y,col,IANd,JANd,PDJ)

  ! full FD column (all reactions)
  dY = eps*max(abs(Y(col)),1.d-30)
  Yp(1:N)=Y(1:N); Yp(col)=Y(col)+dY ; Ym(1:N)=Y(1:N); Ym(col)=Y(col)-dY
  call get_temporal_derivatives(N,0.d0,Yp,ydotp)
  call get_temporal_derivatives(N,0.d0,Ym,ydotm)
  fd(1:N) = (ydotp(1:N)-ydotm(1:N))/(2.d0*dY)

  ! hand-summed type-50 Jacobian wrt Y(col): for every ice-transport reaction,
  ! d(flux)/dY(col) is nonzero only where col is a reactant (slot 1 or 2).
  ana_t50(1:N) = 0.d0
  do rr=1,nb_reactions
    if (REACTION_TYPE(rr).ne.COAGULATION_TYPE) cycle
    if (REACTION_COMPOUNDS_ID(1,rr).lt.1 .or. REACTION_COMPOUNDS_ID(1,rr).gt.N) cycle
    if (SPECIES_PHASE(REACTION_COMPOUNDS_ID(1,rr)).ne.1) cycle       ! ice-transport only
    df = 0.d0
    if (REACTION_COMPOUNDS_ID(1,rr).eq.col) df = df + reaction_rates(rr)*Y(REACTION_COMPOUNDS_ID(2,rr))*actual_gas_density
    if (REACTION_COMPOUNDS_ID(2,rr).eq.col) df = df + reaction_rates(rr)*Y(REACTION_COMPOUNDS_ID(1,rr))*actual_gas_density
    if (df.eq.0.d0) cycle
    do s=4,8
      cid=REACTION_COMPOUNDS_ID(s,rr)
      if (cid.ge.1.and.cid.le.N) ana_t50(cid)=ana_t50(cid)+REACTION_PRODUCT_WEIGHTS(s,rr)*df
    enddo
    do s=1,3
      cid=REACTION_COMPOUNDS_ID(s,rr)
      if (cid.ge.1.and.cid.le.N) ana_t50(cid)=ana_t50(cid)-df
    enddo
  enddo

  write(*,'(/,a)') '=== STATE C (all ice 1e-10), column GRAIN01^0 ==='
  write(*,'(a)') 'rows where full get_jacobian disagrees with full FD:'
  write(*,'(a)') '  species        hand-t50         get_jac          FD          | t50-vs-FD'
  do row=1,N
    if (abs(PDJ(row)).lt.1.d-28 .and. abs(fd(row)).lt.1.d-28) cycle
    relerr = abs(PDJ(row)-fd(row))/max(abs(fd(row)),1.d-30)
    if (relerr.gt.1.d-3) then
      write(*,'(2x,a12,3es16.6,a,es10.2)') species_name(row), ana_t50(row), PDJ(row), fd(row), &
        '  |', abs(ana_t50(row)-fd(row))/max(abs(fd(row)),1.d-30)
    endif
  enddo

  ! decompose one worst row: list every ice-transport reaction contributing to it
  write(*,'(/,a)') '=== decomposition of row J03CH2OH in column GRAIN01^0 ==='
  do row=1,N
    if (trim(species_name(row)).ne.'J03CH2OH') cycle
    do rr=1,nb_reactions
      if (REACTION_TYPE(rr).ne.COAGULATION_TYPE) cycle
      if (REACTION_COMPOUNDS_ID(1,rr).lt.1 .or. REACTION_COMPOUNDS_ID(1,rr).gt.N) cycle
      if (SPECIES_PHASE(REACTION_COMPOUNDS_ID(1,rr)).ne.1) cycle
      df = 0.d0
      if (REACTION_COMPOUNDS_ID(1,rr).eq.col) df = df + reaction_rates(rr)*Y(REACTION_COMPOUNDS_ID(2,rr))*actual_gas_density
      if (REACTION_COMPOUNDS_ID(2,rr).eq.col) df = df + reaction_rates(rr)*Y(REACTION_COMPOUNDS_ID(1,rr))*actual_gas_density
      if (df.eq.0.d0) cycle
      ! net weight of this row in this reaction
      block
        real(double_precision) :: netw ; integer :: s
        netw=0.d0
        do s=4,8
          if (REACTION_COMPOUNDS_ID(s,rr).eq.row) netw=netw+REACTION_PRODUCT_WEIGHTS(s,rr)
        enddo
        do s=1,3
          if (REACTION_COMPOUNDS_ID(s,rr).eq.row) netw=netw-1.d0
        enddo
        if (netw.ne.0.d0) then
          call print_reaction(rr)
          write(*,'(a,es14.6,a,es14.6)') '        contributes netw*df = ', netw*df, '   (netw=',netw
        endif
      end block
    enddo
  enddo
  deallocate(ana_t50)
end block

! ======================================================================
! (10) Localize the FD source for one row: compare get_temporal_derivatives FD
!      against a full hand-reconstructed RHS (sum over ALL reactions of
!      flux = rate*prod(Y_reactants)*density^(nbody-1), deposited by weights).
!      If the two FDs agree, the extra coupling is a real reaction-flux term my
!      per-reaction t50 scan mis-enumerated; if they disagree, get_temporal_
!      derivatives injects a coupling that is not a plain mass-action flux.
! ======================================================================
block
  integer :: trow, rr, s, c1, c2, c3
  real(double_precision) :: fdlib, fdhand, rr_flux, densq
  real(double_precision), allocatable :: hp(:), hm(:)
  allocate(hp(N), hm(N))

  reaction_rates(1:nb_reactions) = rates_save(1:nb_reactions)
  Y(1:N) = abundances(1:N)
  do row=1,N
    if (SPECIES_PHASE(row).eq.1 .and. Y(row).lt.1.d-10) Y(row)=1.d-10
  enddo
  trow = 0
  do row=1,N
    if (trim(species_name(row)).eq.'J03CH2OH') trow=row
  enddo

  dY = eps*max(abs(Y(col)),1.d-30)
  Yp(1:N)=Y(1:N); Yp(col)=Y(col)+dY ; Ym(1:N)=Y(1:N); Ym(col)=Y(col)-dY

  ! library RHS FD
  call get_temporal_derivatives(N,0.d0,Yp,ydotp)
  call get_temporal_derivatives(N,0.d0,Ym,ydotm)
  fdlib = (ydotp(trow)-ydotm(trow))/(2.d0*dY)

  ! hand RHS FD (all reactions, plain mass action, using the SAME reaction_rates)
  reaction_rates(1:nb_reactions) = rates_save(1:nb_reactions)
  call hand_rhs(Yp, hp) ; call hand_rhs(Ym, hm)
  fdhand = (hp(trow)-hm(trow))/(2.d0*dY)

  write(*,'(/,a)') '=== (10) J03CH2OH FD: library RHS vs hand-reconstructed RHS ==='
  write(*,'(a,es16.6)') '  get_temporal_derivatives FD = ', fdlib
  write(*,'(a,es16.6)') '  hand mass-action RHS FD     = ', fdhand
  write(*,'(a,es16.6)') '  get_jacobian analytic       = ', PDJ(trow)

  ! Which reactions' flux actually changes when Y(col) is perturbed? (hand fluxes)
  write(*,'(a)') '  reactions whose hand-flux changes with Y(col) and touch J03CH2OH:'
  densq = actual_gas_density
  do rr=1,nb_reactions
    c1=REACTION_COMPOUNDS_ID(1,rr); c2=REACTION_COMPOUNDS_ID(2,rr); c3=REACTION_COMPOUNDS_ID(3,rr)
    ! does it touch trow?
    s=0
    if (c1.eq.trow.or.c2.eq.trow.or.c3.eq.trow) s=1
    do row=4,8
      if (REACTION_COMPOUNDS_ID(row,rr).eq.trow) s=1
    enddo
    if (s.eq.0) cycle
    ! does its flux depend on Y(col)?
    if (c1.ne.col.and.c2.ne.col.and.c3.ne.col) cycle
    write(*,'(a,i6,a,i3,a,es11.3)') '    r=',rr,' type=',REACTION_TYPE(rr),' rate=',reaction_rates(rr)
  enddo

  ! EXACT per-reaction attribution of Delta(dY_trow/dt) between Yp and Ym.
  write(*,'(a)') '  per-reaction contribution to the FD of J03CH2OH (nonzero only):'
  do rr=1,nb_reactions
    block
      real(double_precision) :: fp, fm, netw, contrib
      integer :: ss, cc
      c1=REACTION_COMPOUNDS_ID(1,rr); c2=REACTION_COMPOUNDS_ID(2,rr); c3=REACTION_COMPOUNDS_ID(3,rr)
      if (c1.lt.1.or.c1.gt.N) cycle
      ! flux at Yp / Ym
      if (c2.ge.1.and.c2.le.N) then
        if (c3.ge.1.and.c3.le.N) then
          fp=reaction_rates(rr)*Yp(c1)*Yp(c2)*Yp(c3)*actual_gas_density*actual_gas_density
          fm=reaction_rates(rr)*Ym(c1)*Ym(c2)*Ym(c3)*actual_gas_density*actual_gas_density
        else
          fp=reaction_rates(rr)*Yp(c1)*Yp(c2)*actual_gas_density
          fm=reaction_rates(rr)*Ym(c1)*Ym(c2)*actual_gas_density
        endif
      else
        fp=reaction_rates(rr)*Yp(c1) ; fm=reaction_rates(rr)*Ym(c1)
      endif
      netw=0.d0
      do ss=4,8
        if (REACTION_COMPOUNDS_ID(ss,rr).eq.trow) netw=netw+REACTION_PRODUCT_WEIGHTS(ss,rr)
      enddo
      do ss=1,3
        if (REACTION_COMPOUNDS_ID(ss,rr).eq.trow) netw=netw-1.d0
      enddo
      contrib=netw*(fp-fm)/(2.d0*dY)
      if (abs(contrib).gt.1.d-30) then
        call print_reaction(rr)
        write(*,'(a,es14.6,a,f7.3,a)') '        FD-contrib=',contrib,'  (netw=',netw,')'
      endif
    end block
  enddo
  deallocate(hp,hm)
end block

! ======================================================================
! (11) Confirm the FD failure is background cancellation, not a Jacobian bug:
!      (a) show the magnitude of dY_J03CH2OH/dt (the Y(GRAIN01)-independent
!          chemistry background that the subtraction must cancel), and
!      (b) recompute the FD with all NON-type-50 rates zeroed (ice-transport
!          only, no large background) -- it should then equal the analytic.
! ======================================================================
block
  integer :: trow, rr
  real(double_precision) :: fd_coag
  real(double_precision), allocatable :: hp(:), hm(:), rz(:)
  allocate(hp(N), hm(N), rz(nb_reactions))
  Y(1:N) = abundances(1:N)
  do row=1,N
    if (SPECIES_PHASE(row).eq.1 .and. Y(row).lt.1.d-10) Y(row)=1.d-10
  enddo
  trow=0 ; do row=1,N ; if (trim(species_name(row)).eq.'J03CH2OH') trow=row ; enddo
  dY = eps*max(abs(Y(col)),1.d-30)
  Yp(1:N)=Y(1:N); Yp(col)=Y(col)+dY ; Ym(1:N)=Y(1:N); Ym(col)=Y(col)-dY

  reaction_rates(1:nb_reactions) = rates_save(1:nb_reactions)
  call hand_rhs(Yp, hp) ; call hand_rhs(Ym, hm)
  write(*,'(/,a)') '=== (11) cancellation check for J03CH2OH ==='
  write(*,'(a,2es16.6)') '  dY_J03CH2OH/dt at Yp, Ym (background)      = ', hp(trow), hm(trow)
  write(*,'(a,es16.6)')  '  their difference hp-hm                     = ', hp(trow)-hm(trow)
  write(*,'(a,es16.6)')  '  true signal (analytic)*2dY                 = ', 5.53d-17*2.d0*dY
  write(*,'(a,es10.2)')  '  => signal/background ratio ~ ', abs(5.53d-17)/max(abs(hp(trow)),1.d-30)

  ! ice-transport only: zero every non-type-50 rate, keep the constant kernel
  rz(1:nb_reactions) = rates_save(1:nb_reactions)
  do rr=1,nb_reactions
    if (REACTION_TYPE(rr).ne.COAGULATION_TYPE) rz(rr)=0.d0
  enddo
  reaction_rates(1:nb_reactions) = rz(1:nb_reactions)
  call hand_rhs(Yp, hp) ; call hand_rhs(Ym, hm)
  fd_coag = (hp(trow)-hm(trow))/(2.d0*dY)
  write(*,'(a,es16.6)') '  FD with chemistry zeroed (ice-transport only)= ', fd_coag
  write(*,'(a,es16.6)') '  analytic get_jacobian                        = ', 5.53d-17
  reaction_rates(1:nb_reactions) = rates_save(1:nb_reactions)
  deallocate(hp,hm,rz)
end block

! ======================================================================
! (12) Same treatment for a grain-grain coag row (GRAIN02) seen with
!      get_jac != 0 but full-FD = 0, to classify it as cancellation or real.
! ======================================================================
block
  integer :: trow, rr, ss
  real(double_precision) :: fd_full, fd_attr, contrib, fp, fm
  real(double_precision), allocatable :: hp(:), hm(:)
  integer :: c1,c2,c3
  allocate(hp(N),hm(N))
  reaction_rates(1:nb_reactions) = rates_save(1:nb_reactions)
  Y(1:N) = abundances(1:N)
  do row=1,N ; if (SPECIES_PHASE(row).eq.1 .and. Y(row).lt.1.d-10) Y(row)=1.d-10 ; enddo
  trow=0 ; do row=1,N ; if (trim(species_name(row)).eq.'GRAIN02') trow=row ; enddo
  dY = eps*max(abs(Y(col)),1.d-30)
  Yp(1:N)=Y(1:N); Yp(col)=Y(col)+dY ; Ym(1:N)=Y(1:N); Ym(col)=Y(col)-dY
  call hand_rhs(Yp,hp) ; call hand_rhs(Ym,hm)
  fd_full = (hp(trow)-hm(trow))/(2.d0*dY)
  fd_attr = 0.d0
  do rr=1,nb_reactions
    c1=REACTION_COMPOUNDS_ID(1,rr); c2=REACTION_COMPOUNDS_ID(2,rr); c3=REACTION_COMPOUNDS_ID(3,rr)
    if (c1.lt.1.or.c1.gt.N) cycle
    if (c2.ge.1.and.c2.le.N) then
      if (c3.ge.1.and.c3.le.N) then
        fp=reaction_rates(rr)*Yp(c1)*Yp(c2)*Yp(c3)*actual_gas_density*actual_gas_density
        fm=reaction_rates(rr)*Ym(c1)*Ym(c2)*Ym(c3)*actual_gas_density*actual_gas_density
      else
        fp=reaction_rates(rr)*Yp(c1)*Yp(c2)*actual_gas_density
        fm=reaction_rates(rr)*Ym(c1)*Ym(c2)*actual_gas_density
      endif
    else
      fp=reaction_rates(rr)*Yp(c1); fm=reaction_rates(rr)*Ym(c1)
    endif
    contrib=0.d0
    do ss=4,8 ; if (REACTION_COMPOUNDS_ID(ss,rr).eq.trow) contrib=contrib+REACTION_PRODUCT_WEIGHTS(ss,rr) ; enddo
    do ss=1,3 ; if (REACTION_COMPOUNDS_ID(ss,rr).eq.trow) contrib=contrib-1.d0 ; enddo
    fd_attr = fd_attr + contrib*(fp-fm)/(2.d0*dY)
  enddo
  write(*,'(/,a)') '=== (12) GRAIN02 row, column GRAIN01^0 ==='
  write(*,'(a,2es16.6)') '  dY_GRAIN02/dt at Yp,Ym (background)   = ', hp(trow), hm(trow)
  write(*,'(a,es16.6)')  '  full hand FD                          = ', fd_full
  write(*,'(a,es16.6)')  '  per-reaction FD attribution (no accum)= ', fd_attr
  deallocate(hp,hm)
end block

call exit(0)

contains

  subroutine print_reaction(rr)
  integer, intent(in) :: rr
  integer :: ci, s
  character(len=160) :: line
  character(len=16) :: tok
  line = ''
  do s=1,3
    ci = REACTION_COMPOUNDS_ID(s,rr)
    if (ci.ge.1 .and. ci.le.N) then
      if (len_trim(line).gt.0) line = trim(line)//' + '
      line = trim(line)//trim(species_name(ci))
    endif
  enddo
  line = trim(line)//'  ->  '
  do s=4,8
    ci = REACTION_COMPOUNDS_ID(s,rr)
    if (ci.ge.1 .and. ci.le.N) then
      write(tok,'(f5.3,1x)') REACTION_PRODUCT_WEIGHTS(s,rr)
      line = trim(line)//' '//trim(adjustl(tok))//trim(species_name(ci))
    endif
  enddo
  write(*,'(a,i6,a,i3,a,es11.3,a,a)') ' r=',rr,' type=',REACTION_TYPE(rr), &
      ' rate=',reaction_rates(rr),'  ',trim(line)
  end subroutine print_reaction

  !> Single-reaction RHS: flux = rate*Y(react1)*Y(react2)*density (two-body),
  !! deposited +w_p to products and -1 to reactants. Mirrors get_temporal_derivatives
  !! for one two-body reaction.
  subroutine rhs_single(rr, rate, cc, ww, Yin, out)
  integer, intent(in) :: rr, cc(8)
  real(double_precision), intent(in) :: rate, ww(8), Yin(:)
  real(double_precision), intent(out) :: out(:)
  real(double_precision) :: flux
  integer :: s
  out(1:N) = 0.d0
  if (cc(2).ge.1 .and. cc(2).le.N) then
    flux = rate*Yin(cc(1))*Yin(cc(2))*actual_gas_density   ! two-body (ice transport)
  else
    flux = rate*Yin(cc(1))                        ! one-body (not expected for ice transport)
  endif
  do s=4,8
    if (cc(s).ge.1 .and. cc(s).le.N) out(cc(s)) = out(cc(s)) + ww(s)*flux
  enddo
  do s=1,3
    if (cc(s).ge.1 .and. cc(s).le.N) out(cc(s)) = out(cc(s)) - flux
  enddo
  end subroutine rhs_single

  !> Full RHS reconstructed as plain mass action over ALL reactions, using the
  !! current reaction_rates: flux = rate * prod(Y_reactants) * density^(nbody-1),
  !! deposited +w_p to products, -flux to each reactant. No SUMLAY/GTODN re-read.
  subroutine hand_rhs(Yin, out)
  real(double_precision), intent(in) :: Yin(:)
  real(double_precision), intent(out) :: out(:)
  real(double_precision) :: flux
  integer :: rr, s, c1, c2, c3, cc
  out(1:N) = 0.d0
  do rr=1,nb_reactions
    c1=REACTION_COMPOUNDS_ID(1,rr); c2=REACTION_COMPOUNDS_ID(2,rr); c3=REACTION_COMPOUNDS_ID(3,rr)
    if (c1.lt.1 .or. c1.gt.N) cycle
    if (c2.ge.1 .and. c2.le.N) then
      if (c3.ge.1 .and. c3.le.N) then
        flux = reaction_rates(rr)*Yin(c1)*Yin(c2)*Yin(c3)*actual_gas_density*actual_gas_density
      else
        flux = reaction_rates(rr)*Yin(c1)*Yin(c2)*actual_gas_density
      endif
    else
      flux = reaction_rates(rr)*Yin(c1)
    endif
    do s=4,8
      cc=REACTION_COMPOUNDS_ID(s,rr)
      if (cc.ge.1 .and. cc.le.N) out(cc)=out(cc)+REACTION_PRODUCT_WEIGHTS(s,rr)*flux
    enddo
    do s=1,3
      cc=REACTION_COMPOUNDS_ID(s,rr)
      if (cc.ge.1 .and. cc.le.N) out(cc)=out(cc)-flux
    enddo
  enddo
  end subroutine hand_rhs

end program ice_transport_jac_diag
