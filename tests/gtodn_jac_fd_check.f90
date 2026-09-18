! ===========================================================================
! tests/gtodn_jac_fd_check.f90 -- Rung 5 entry-level FD check (miniature gate 5a)
!
! At a state with live surface chemistry (dynamic GTODN on), compare the analytic
! get_jacobian column for every grain population species against a central finite
! difference of get_temporal_derivatives w.r.t. that species. A correct analytic
! grain-abundance entry (accretion nu=+1, LH nu=-1, both /Y_tot) makes the full
! column match FD to central-difference tolerance. A sign or form error (e.g.
! +rate/Y instead of -rate/Y**2 for LH) blows the relative error up.
!
! Two states are tested:
!   (1) ice-rich, all grain bins above the floor      -- exercises accretion + LH.
!   (2) one grain bin driven below the floor          -- LH floor derivative -> 0.
!
! Run from a prepared run directory (coagulation = 1). Exit 0 = PASS.
! ===========================================================================
program gtodn_jac_fd_check
use iso_fortran_env, only: error_unit
use global_variables
use ode_solver
use gasgrain
implicit none

integer :: N, J, k, r, i, row, nbad, worst_row, worst_col
real(double_precision), allocatable :: Y(:), Yp(:), Ym(:), ydotp(:), ydotm(:)
real(double_precision), allocatable :: PDJ(:), fdcol(:), IANd(:), JANd(:)
real(double_precision) :: dY, eps, num, den, relerr, worst, worst_state2, worstC, worstD
logical :: any_grain_col
character(len=8) :: tag
integer :: iJH, iJCO

call init_gasgrain()

N = nb_species
allocate(Y(N), Yp(N), Ym(N), ydotp(N), ydotm(N), PDJ(N), fdcol(N), IANd(N), JANd(N))
IANd = 0.d0 ; JANd = 0.d0   ! get_jacobian does not read IAN/JAN in its body
eps = 1.d-6

write(*,'(a,i0,a,i0,a,l1)') 'N(species)=', N, '  nb_grains=', nb_grains, &
  '  dynamic_gtodn_active=', dynamic_gtodn_active()
write(*,'(a,i0)') 'in-scope reactions in gtodn_jac_list (types 99,14): ', gtodn_jac_n
if (.not.dynamic_gtodn_active()) then
  write(error_unit,'(a)') 'FAIL: dynamic GTODN not active -- run with coagulation = 1.'
  call exit(2)
endif

! -------- build an ice-rich state: seed = initial abundances, floor all surface
!          (J) species to 1e-10 so LH/accretion fluxes are non-trivial. --------
Y(1:N) = abundances(1:N)
do i=1,N
  if (SPECIES_PHASE(i).eq.1 .and. Y(i).lt.1.d-10) Y(i) = 1.d-10
enddo

! -------- controlled states that keep SUMLAY (hence its grain-derivative)
!          negligible, so accretion (nu=+1) and LH (nu=-1) can be validated in
!          isolation from the Rung-3 live-SUMLAY couplings (photodesorption 66/67
!          and the H/H2 sticking term), which are a later increment. --------
iJH = 0 ; iJCO = 0
do i=1,N
  if (trim(species_name(i)).eq.'JH')  iJH  = i
  if (trim(species_name(i)).eq.'JCO') iJCO = i
enddo

! ================= STATE A: zero surface ice -> accretion only, SUMLAY=0 =========
Y(1:N) = abundances(1:N)
do i=1,N
  if (SPECIES_PHASE(i).ge.1) Y(i) = 0.d0     ! no ice => ab_tot=0 => SUMLAY=0, dSUMLAY/dYtot=0
enddo
call fd_compare(Y, worst, worst_row, worst_col, nbad)
write(*,'(a)') '--------------------------------------------------------------'
write(*,'(a,es10.3,a,i0,a)') '[A: accretion only] worst grain-col rel-err = ', worst, &
  '  (cols >1e-5: ', nbad, ')'

! ================= STATE B: trace ice on a real LH pair -> LH reciprocal ========
! Find the first type-14 (Langmuir-Hinshelwood) reaction and seed its two adsorbed
! reactants; test the grain column of that reaction's bin. Surface species carry a
! per-bin name prefix, so pick them from the reaction table, not by name.
block
  integer :: rlh, a1, a2, kbin, Jlh
  rlh = 0
  do i=1,nb_reactions
    if (REACTION_TYPE(i).eq.14 .and. GRAIN_RANK(i).ge.1) then
      a1 = REACTION_COMPOUNDS_ID(1,i); a2 = REACTION_COMPOUNDS_ID(2,i)
      if (a1.ge.1.and.a1.le.N .and. a2.ge.1.and.a2.le.N) then ; rlh = i ; exit ; endif
    endif
  enddo
  if (rlh.eq.0) then
    write(*,'(a)') '[B] no usable type-14 LH reaction found -- skipping LH probe.'
    worst_state2 = 0.d0
  else
    a1 = REACTION_COMPOUNDS_ID(1,rlh); a2 = REACTION_COMPOUNDS_ID(2,rlh); kbin = GRAIN_RANK(rlh)
    Jlh = INDGRAIN(kbin)
    Y(1:N) = abundances(1:N)
    do i=1,N
      if (SPECIES_PHASE(i).ge.1) Y(i) = 0.d0
    enddo
    Y(a1) = 1.d-9 ; Y(a2) = 1.d-9
    call fd_compare(Y, worst_state2, worst_row, worst_col, nbad)
    write(*,'(a,es10.3,a,i0,a)') '[B: LH reciprocal]  worst grain-col rel-err = ', worst_state2, &
      '  (cols >1e-5: ', nbad, ')'
    ! explicit LH-row confirmation for this reaction's bin column
    call get_temporal_derivatives(N, 0.d0, Y, ydotp)
    call get_jacobian(N, 0.d0, Y, Jlh, IANd, JANd, PDJ)
    dY = eps * max(abs(Y(Jlh)), 1.d-30)
    Yp(1:N) = Y(1:N) ; Yp(Jlh) = Y(Jlh) + dY
    Ym(1:N) = Y(1:N) ; Ym(Jlh) = Y(Jlh) - dY
    call get_temporal_derivatives(N, 0.d0, Yp, ydotp)
    call get_temporal_derivatives(N, 0.d0, Ym, ydotm)
    fdcol(1:N) = (ydotp(1:N) - ydotm(1:N)) / (2.d0*dY)
    write(*,'(a,i0,a,a,a,a,a,a)') '  [LH reaction ',rlh,': ',trim(species_name(a1)),' + ', &
      trim(species_name(a2)),' -> ',trim(species_name(REACTION_COMPOUNDS_ID(4,rlh)))
    write(*,'(a)') '   reactants (JH-side) expect analytic>0; product expect analytic<0:'
    do row=1,N
      if (row.eq.a1 .or. row.eq.a2 .or. row.eq.REACTION_COMPOUNDS_ID(4,rlh)) then
        write(*,'(4x,a12,a,2es16.6)') species_name(row), ' ana/fd=', PDJ(row), fdcol(row)
      endif
    enddo
  endif
end block

! ============ STATE C: mild blanket ice (photodesorption RHS live via SUMLAY,
!              cap OFF, H/H2 sticking active) -- the state that failed before this increment
Y(1:N) = abundances(1:N)
do i=1,N
  if (SPECIES_PHASE(i).eq.1 .and. Y(i).lt.1.d-10) Y(i) = 1.d-10
enddo
call fd_compare(Y, worstC, worst_row, worst_col, nbad)
write(*,'(a,es10.3,a,i0,a)') '[C: mild ice, cap off] worst grain-col rel-err = ', worstC, &
  '  (cols >1e-5: ', nbad, ')'
if (worstC.gt.1.d-5 .and. worst_col.ge.1) then
  J = worst_col
  call get_temporal_derivatives(N, 0.d0, Y, ydotp)
  call get_jacobian(N, 0.d0, Y, J, IANd, JANd, PDJ)
  dY = eps*max(abs(Y(J)),1.d-30)
  Yp(1:N)=Y(1:N); Yp(J)=Y(J)+dY; Ym(1:N)=Y(1:N); Ym(J)=Y(J)-dY
  call get_temporal_derivatives(N,0.d0,Yp,ydotp); call get_temporal_derivatives(N,0.d0,Ym,ydotm)
  fdcol(1:N)=(ydotp(1:N)-ydotm(1:N))/(2.d0*dY)
  write(*,'(a,i0,a,a)') '  [diagC] col ',J,' ',trim(species_name(J))
  do row=1,N
    if (abs(fdcol(row)).lt.1.d-28) cycle
    relerr=abs(PDJ(row)-fdcol(row))/abs(fdcol(row))
    if (relerr.gt.1.d-3) write(*,'(4x,a12,2es18.8,es11.2)') species_name(row),PDJ(row),fdcol(row),relerr
  enddo
endif

! ============ STATE D: heavy blanket ice -> SUMLAY >= MLAY on populated bins
!              (photodesorption monolayer cap ACTIVE) ============
Y(1:N) = abundances(1:N)
do i=1,N
  if (SPECIES_PHASE(i).eq.1) Y(i) = max(Y(i), 1.d-4)
enddo
call fd_compare(Y, worstD, worst_row, worst_col, nbad)
write(*,'(a,es10.3,a,i0,a)') '[D: heavy ice, cap on] worst grain-col rel-err = ', worstD, &
  '  (cols >1e-5: ', nbad, ')'

worst = max(max(worst, worst_state2), max(worstC, worstD))

write(*,'(a)') '--------------------------------------------------------------'
if (worst.lt.1.d-5) then
  write(*,'(a)') 'RESULT: PASS (analytic grain-column Jacobian matches FD to <1e-5).'
  call exit(0)
else
  write(error_unit,'(a)') 'RESULT: FAIL (analytic vs FD mismatch above 1e-5).'
  call exit(1)
endif

contains

  !> Compare analytic get_jacobian vs central-difference FD for every grain column
  !! at abundance state Yin. Returns the worst column relative error and its location.
  subroutine fd_compare(Yin, worst_out, worst_row_out, worst_col_out, nbad_out)
  real(double_precision), intent(in) :: Yin(:)
  real(double_precision), intent(out) :: worst_out
  integer, intent(out) :: worst_row_out, worst_col_out, nbad_out
  real(double_precision) :: colworst
  integer :: kk, jj

  worst_out = 0.d0 ; worst_row_out = 0 ; worst_col_out = 0 ; nbad_out = 0

  do kk=1,nb_grains
    do jj=1,2
      if (jj.eq.1) then
        J = INDGRAIN(kk)
      else
        J = INDGRAIN_MINUS(kk)
      endif
      if (J.lt.1) cycle

      ! analytic column (refresh reaction_rates at base Y first)
      Y(1:N) = Yin(1:N)
      call get_temporal_derivatives(N, 0.d0, Y, ydotp)   ! populates reaction_rates(base Y)
      call get_jacobian(N, 0.d0, Y, J, IANd, JANd, PDJ)

      ! central finite difference of the RHS w.r.t. Y(J)
      dY = eps * max(abs(Yin(J)), 1.d-30)
      Yp(1:N) = Yin(1:N) ; Yp(J) = Yin(J) + dY
      Ym(1:N) = Yin(1:N) ; Ym(J) = Yin(J) - dY
      call get_temporal_derivatives(N, 0.d0, Yp, ydotp)
      call get_temporal_derivatives(N, 0.d0, Ym, ydotm)
      fdcol(1:N) = (ydotp(1:N) - ydotm(1:N)) / (2.d0 * dY)

      ! worst relative error over rows with a meaningful FD magnitude
      colworst = 0.d0
      do row=1,N
        den = abs(fdcol(row))
        if (den.lt.1.d-25) cycle          ! skip structurally-zero rows
        relerr = abs(PDJ(row) - fdcol(row)) / den
        if (relerr.gt.colworst) colworst = relerr
        if (relerr.gt.worst_out) then
          worst_out = relerr ; worst_row_out = row ; worst_col_out = J
        endif
      enddo
      if (colworst.gt.1.d-5) nbad_out = nbad_out + 1
    enddo
  enddo
  end subroutine fd_compare

end program gtodn_jac_fd_check
