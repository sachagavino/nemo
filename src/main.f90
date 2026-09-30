! -----------------------------------------------------------------------
!  0D multi-grain gas-grain chemical code.
!
!  Skeleton extracted from NMGC-2.0 (github.com/sachagavino/nmgc-2.0),
!  itself derived from Nautilus (Hersant, Wakelam, Cossou, Gavino, Iqbal).
!
!  What this driver does, and only this:
!    - initialise the model (init_gasgrain)
!    - loop on output times
!    - integrate the chemical scheme with DLSODES between two output times
!    - write the outputs
!
!  Everything 1D (spatial grid, species diffusion, Crank-Nicholson,
!  1D_static.dat) has been removed, as have the on-the-fly dust temperature
!  computation, the disk photorates, and the rates/major-reactions
!  post-processing commands.
!
!  HOOK (coagulation): the dust size distribution will be integrated INSIDE
!  the ODE system, not operator-split around it. There is therefore no
!  sub-stepping loop here on purpose: adding coagulation must not require
!  touching this driver, only the RHS/Jacobian and the state vector.
!
!  INPUT FILES
!    parameters.in           : parameters and flags of the model
!    abundances.in           : initial abundances (gas + ice)
!    dust_grid_table.in      : tabulated grain bins (radius, 1/abundance, Td, T_CR,peak);
!                              only when dust_grid_source/dust_ic = tabulated
!                              -> disappears at stage 3, replaced by a derived
!                                 grid + a separate dust IC file
!    element.in              : name and mass [AMU] of the elements
!    gas_species.in / gas_reactions.in
!    grain_species.in / grain_reactions.in
!    activation_energies.in  : activation energies of endothermic reactions
!    surface_parameters.in   : energies and parameters of surface processes
!
!  OUTPUT FILES
!    abundances.out          : binary, one record set per output time
!    rates.out               : binary, reaction rates at each output time
!    species.out             : species <-> index correspondence
!    elemental_abundances.out/.tmp
!    abundances.tmp          : ASCII snapshot of the last output time
! -----------------------------------------------------------------------

PROGRAM nmgc

  use global_variables
  use iso_fortran_env
  use shielding
  use dust_evolution
  use utilities
  use gasgrain
  use outputs

  implicit none

  ! Parameters for DLSODES. RTOL is the RELATIVE_TOLERANCE parameter in global_variables.f90
  integer :: itol = 2 !< ITOL = 1 or 2 according as ATOL is a scalar or array.
  integer :: itask = 1 !< ITASK = 1 for normal computation of output values of Y at t = TOUT.
  integer :: istate = 1 !< ISTATE = integer flag (input and output). Set ISTATE = 1.
  integer :: iopt = 1 !< IOPT = 1 to indicate optional inputs are used.
  integer :: mf !< DLSODES method flag; set by solver_method_flag() after init (coag-on+symbolic -> 021 analytic; coag-off or numerical -> 121; both MITER=1. FD-complete 022 only via NEMO_ORACLE_MF)
  real(double_precision) :: atol = 1.d-99 !< absolute tolerance parameter

  real(double_precision) :: output_timestep !< Timestep to reach the next output time [s]

  real(double_precision), dimension(:), allocatable :: output_times !< [s] dim(NB_OUTPUTS)
  real(double_precision) :: output_step !< [s] linear or log spacing of the output times
  integer :: output_idx !< index of the output time loop

  real(double_precision) :: code_start_time, code_current_time, code_elapsed_time
  real(double_precision) :: remaining_time !< estimated remaining time [s]

  ! Rung 4 solver characterization (4-IV cost + convergence effort): accumulated across all DLSODES calls.
  integer(kind=8) :: tot_nst = 0  !< total accepted steps (IWORK(11))
  integer(kind=8) :: tot_nfe = 0  !< total RHS evaluations (IWORK(12))
  integer(kind=8) :: tot_nje = 0  !< total Jacobian evals / LU decompositions (IWORK(13))
  integer(kind=8) :: tot_fail = 0 !< number of DLSODES calls returning istate /= 2

  integer :: i, ic_i !< For loops
  character(2) :: c_i !< to convert the grain rank encoded in a species name into an integer
  logical :: gotit, quit
  logical :: isDefined

  stdo = 6
  ffli = 5

  do_chemical_scheme = .false.
  do_outputs         = .false.

  !================================================================
  !                          START ACTION
  !================================================================

  call interpet_command_line(gotit,.false.,quit)
  if (.not.gotit) goto 650

  if (do_chemical_scheme) then

    call write_banner()

    call init_gasgrain()

    ! Select the DLSODES method flag now that the coagulation switch is known.
    mf = solver_method_flag()

    call initialize_work_arrays()

    ! ---- list of output times
    select case(OUTPUT_TYPE)
      case('linear')
        allocate(output_times(NB_OUTPUTS))
        output_step = (STOP_TIME - START_TIME) / dfloat(NB_OUTPUTS - 1)
        do i=1, NB_OUTPUTS - 1
          output_times(i) = START_TIME + output_step * dfloat(i)
        enddo
        output_times(NB_OUTPUTS) = STOP_TIME ! To ensure the exact same final value

      case('log')
        allocate(output_times(NB_OUTPUTS))
        output_step = (STOP_TIME/START_TIME) ** (1.d0/dfloat(NB_OUTPUTS - 1))
        do i=1, NB_OUTPUTS - 1
          output_times(i) = START_TIME * output_step ** (i - 1.d0)
        enddo
        output_times(NB_OUTPUTS) = STOP_TIME ! To ensure the exact same final value

      case default
        write(error_unit,*) 'The OUTPUT_TYPE="', trim(OUTPUT_TYPE),'" cannot be found.'
        write(error_unit,*) 'Values possible : linear, log'
        write(error_unit,'(a)') 'Error in main.'
        call exit(12)
    end select

    ! Phase III: full active reaction range (coupled production path). The split driver
    ! narrows it per sub-step and restores it; nothing else in the code writes it.
    active_lo = 1
    active_hi = nb_reactions
    if (split_mode.ne.0) call split_setup_and_audit()

    current_time = 0.d0
    call cpu_time(code_start_time)

    ! ---- loop on output times
    do output_idx=1, NB_OUTPUTS

      output_timestep = output_times(output_idx) - current_time

      call get_structure_properties(time=current_time, & ! Inputs
                              av=visual_extinction, density=H_number_density, & ! Outputs
                              gas_temperature=gas_temperature) ! Outputs

      call get_grain_temperature(time=current_time, gastemperature=gas_temperature, av=visual_extinction, & ! Inputs
                              grain_temperature=dust_temperature) ! Outputs

      ! We can't add input variables in dlsodes-called routines, so we store them as global variables
      actual_gas_temp     = gas_temperature
      actual_dust_temp(:) = dust_temperature(:)
      actual_av           = visual_extinction
      actual_gas_density  = H_number_density

      DO i=1,nb_grains
        IF ((actual_dust_temp(i).gt.15.0d+00).and.(is_er_cir.eq.1)) THEN
          PRINT*, "Warning, for grain",i,"with radius =",grain_radii(i)
          PRINT*, "the temperature range is above the validity domain of the complex induced reaction mechanism..."
        ENDIF
      ENDDO

      ! Column densities of the self-shielded species, from Av. Each is only used
      ! by the corresponding H2/CO/N2 self-shielding rate code, which is absent when
      ! those reactions are not in the network (e.g. the dust-only coagulation
      ! fixture). Guard on a resolved index so a species missing from a reduced
      ! network does not index abundances(0); byte-identical for any full network,
      ! where every index is positive.
      NH=0.d0; NH2=0.d0; NN2=0.d0; NCO=0.d0; NH2O=0.d0; NCH=0.d0; NCH3=0.d0
      NH2CO=0.d0; NCO2=0.d0; NN2O=0.d0; NCH4=0.d0; NOH=0.d0; NHCO=0.d0; NCN=0.d0
      NHCN=0.d0; NHNC=0.d0; NNH=0.d0; NNH2=0.d0; NNH3=0.d0
      if (indH   >0) NH    = actual_av/AV_NH_ratio * abundances(indH)
      if (indH2  >0) NH2   = actual_av/AV_NH_ratio * abundances(indH2)
      if (indN2  >0) NN2   = actual_av/AV_NH_ratio * abundances(indN2)
      if (indCO  >0) NCO   = actual_av/AV_NH_ratio * abundances(indCO)
      if (indH2O >0) NH2O  = actual_av/AV_NH_ratio * abundances(indH2O)
      if (indCH  >0) NCH   = actual_av/AV_NH_ratio * abundances(indCH)
      if (indCH3 >0) NCH3  = actual_av/AV_NH_ratio * abundances(indCH3)
      if (indH2CO>0) NH2CO = actual_av/AV_NH_ratio * abundances(indH2CO)
      if (indCO2 >0) NCO2  = actual_av/AV_NH_ratio * abundances(indCO2)
      if (indN2O >0) NN2O  = actual_av/AV_NH_ratio * abundances(indN2O)
      if (indCH4 >0) NCH4  = actual_av/AV_NH_ratio * abundances(indCH4)
      if (indOH  >0) NOH   = actual_av/AV_NH_ratio * abundances(indOH)
      if (indHCO >0) NHCO  = actual_av/AV_NH_ratio * abundances(indHCO)
      if (indCN  >0) NCN   = actual_av/AV_NH_ratio * abundances(indCN)
      if (indHCN >0) NHCN  = actual_av/AV_NH_ratio * abundances(indHCN)
      if (indHNC >0) NHNC  = actual_av/AV_NH_ratio * abundances(indHNC)
      if (indNH  >0) NNH   = actual_av/AV_NH_ratio * abundances(indNH)
      if (indNH2 >0) NNH2  = actual_av/AV_NH_ratio * abundances(indNH2)
      if (indNH3 >0) NNH3  = actual_av/AV_NH_ratio * abundances(indNH3)

      if (split_mode.eq.0) then
      call integrate_chemical_scheme(delta_t=output_timestep, temp_abundances=abundances(1:nb_species), & ! Inputs
      i_tol=itol, a_tol=atol, i_task=itask, i_opt=iopt, m_f=mf, & ! Inputs
      i_state=istate) ! Output
      else
      ! Phase III operator-split path (isolated; never taken in production).
      call integrate_split(delta_t=output_timestep, temp_abundances=abundances(1:nb_species), &
      i_tol=itol, a_tol=atol, i_task=itask, i_opt=iopt, m_f=mf, i_state=istate)
      endif

      if (istate.eq.-3) stop

      ! Prevent too low abundances, and count the ice layers
      sumlaysurfsave = 0.0d0
      sumlaymantsave = 0.0d0
      do i=1,nb_species
        if (abundances(i).le.1.d-99) then
          abundances(i) = 1.d-99
        endif
        if (species_name(i)(1:1).eq."J") then
          c_i = species_name(i)(2:3)
          read(c_i,'(I2)') ic_i
          sumlaysurfsave(ic_i) = sumlaysurfsave(ic_i) + abundances(i)
        endif
        if (species_name(i)(1:1).eq."K") then
          c_i = species_name(i)(2:3)
          read(c_i,'(I2)') ic_i
          sumlaymantsave(ic_i) = sumlaymantsave(ic_i) + abundances(i)
        endif
      enddo

      do i=1,nb_grains
        sumlaysurfsave(i) = sumlaysurfsave(i) * GTODN(i) / nb_sites_per_grain(i)
        sumlaymantsave(i) = sumlaymantsave(i) * GTODN(i) / nb_sites_per_grain(i)
      enddo

      WRITE(*,'(a)')'----------------------------------------------------------------------------------------------'
      WRITE(*,'(2x,a)') "Physical parameters :"
      WRITE(*,'(4x,a19,ES13.6,1x,a6,4x,a19,ES13.6,1x,a6)') "Density = ",actual_gas_density,"[cm-3]","Av = ",actual_av,"[mag]"
      WRITE(*,'(4x,a16,ES13.6,1x,a6,4x,a27,ES13.6,1x,a10)') "Tgas = ",actual_gas_temp,"[K]", "UV flux = ", UV_flux,"[standard]"
      IF (is_3_phase.eq.0) THEN
        WRITE(*,'(2x,a)') "2 phase model (Rq. Nb. mant. layer should be equal to 0.00):"
      ELSE
        WRITE(*,'(2x,a59,f5.2,a)') "3 phase model (Nb. surf. layer should never be >",nb_active_lay,"):"
      ENDIF

      WRITE(*,'(2x, a)')'Grain parameters :'
      do i=1,nb_grains
        WRITE(*,'(a25,i2, a2)')'Grain rank =',i,': '
        WRITE(*,91) "GTODN = ", GTODN(i),"Grain radius =",grain_radii(i),"[cm]","Grain temperature =",actual_dust_temp(i),"[K]"
        WRITE(*,92) "Nb. surf. layer = ",sumlaysurfsave(i),"Nb. mant. layer = ", sumlaymantsave(i),&
                    &"Nb. total layer = ",sumlaysurfsave(i)+sumlaymantsave(i)
        WRITE(*,'(a)') " "
      enddo
      WRITE(*,'(a)')'=============================================================================================='
  91  FORMAT(2x,a19,x,ES10.3,4x,a14,ES10.2,x,a4,4x,a19,F8.2,x,a4)
  92  FORMAT(2x,a29,x,F10.6,2x,a18,x,F10.6,2x,a18,x,F10.6)

      call check_conservation(abundances(1:nb_species))

      current_time = output_times(output_idx) ! New current time at which abundances are valid

      call cpu_time(code_current_time)
      code_elapsed_time = code_current_time - code_start_time
      remaining_time = code_elapsed_time*NB_OUTPUTS/output_idx - code_elapsed_time

      if (remaining_time.lt.60.d0) then
        write(Output_Unit,'(a,en11.2e2,a,f5.1,a,f5.1,a)') 'T =',current_time/YEAR,&
        ' years [', 100.d0 * output_idx/NB_OUTPUTS, ' %] ; Estimated time remaining: ', remaining_time, ' s'
      else if (remaining_time.lt.3600.d0) then
        write(Output_Unit,'(a,en11.2e2,a,f5.1,a,f5.1,a)') 'T =',current_time/YEAR,&
        ' years [', 100.d0 * output_idx/NB_OUTPUTS, ' %] ; Estimated time remaining: ', remaining_time / 60.d0, ' min.'
      else
        write(Output_Unit,'(a,en11.2e2,a,f5.1,a,f5.1,a)') 'T =',current_time/YEAR,&
        ' years [', 100.d0 * output_idx/NB_OUTPUTS, ' %] ; Estimated time remaining: ', remaining_time / 3600.d0, ' h'
      endif

      call write_current_rates(index=output_idx)
      call write_current_output(index=output_idx)
      ! Overflow guard: the coagulation collision rate the Rung 1 top-bin policy
      ! drops (and the total), recorded per output by write_current_dust. Must stay
      ! ~0; a growing fraction means mass reached the top bins.
      if (coagulation) call dust_coagulation_flux_diag(coag_dropped_flux, coag_total_flux)
      call write_current_dust(index=output_idx)

      first_step_done = .true.
    enddo

    call write_abundances('abundances.tmp')

    ! Rung 4 solver characterization (4-IV cost + convergence effort; oracle vs production).
    write(stdo,'(a)')    ' --- DLSODES solver statistics (Rung 4 characterization) ---'
    write(stdo,'(a,i0)') '   method flag mf              = ', mf
    write(stdo,'(a,i0)') '   steps        NST            = ', tot_nst
    write(stdo,'(a,i0)') '   RHS evals    NFE            = ', tot_nfe
    write(stdo,'(a,i0)') '   Jac/LU evals NJE            = ', tot_nje
    write(stdo,'(a,i0)') '   solver restarts (istate/=2) = ', tot_fail
    write(stdo,'(a,i0)') '   sparse LU    NLU            = ', tot_nlu
    write(stdo,'(a,i0)') '   cold starts  (ISTATE=1)     = ', tot_coldstart
    if (split_mode.ne.0) write(stdo,'(a,i0)') '   Strang macro-steps          = ', tot_macro
    call cpu_time(code_current_time)
    write(stdo,'(a,es12.5)') '   CPU time [s]                = ', code_current_time - code_start_time

  elseif (do_outputs) then
    inquire(file='abundances.out', exist=isDefined)
    if (.not.isDefined) then
      write(Error_unit,'(a)') 'ERROR: there is no output file in the model. Please, compute the chemistry first.'
      stop
    endif

    call write_banner()
    call init_gasgrain()

    write(*,'(a,i0)') 'Number of time outputs: ', nb_outputs
    write(*,'(a,i0)') 'Number of species: ', nb_species
    write(*,'(a,i0)') 'Start generating output ASCII files... '
    call get_outputs()

  end if

  write(stdo,*) 'Done...'
  call flush(stdo)

  650 continue

  contains

! ======================================================================
!> @brief Chemically evolve for a given time delta_t
! ======================================================================
  subroutine integrate_chemical_scheme(delta_t,temp_abundances,i_tol,a_tol,i_task,i_opt,m_f,i_state)
    use global_variables

    implicit none
    ! Inputs
    real(double_precision), intent(in) :: delta_t !<[in] time during which we must integrate
    integer, intent(in) :: i_tol !<[in] ITOL = 1 or 2 according as ATOL is a scalar or array.
    integer, intent(in) :: i_task !<[in] ITASK = 1 for normal computation of output values of Y at t = TOUT.
    integer, intent(in) :: i_opt !<[in] IOPT = 1 to indicate optional inputs are used.
    integer, intent(in) :: m_f !<[in] method flag
    real(double_precision), intent(in) :: a_tol !<[in] integrator tolerance

    ! Outputs
    integer, intent(out) :: i_state !<[out] ISTATE = 2 if DLSODES was successful, negative otherwise.

    ! Input/Output
    real(double_precision), dimension(nb_species), intent(inout) :: temp_abundances !<[in,out] abundances buffer

    ! Locals
    real(double_precision) :: t !< The local time, starting from 0 to delta_t
    real(double_precision), dimension(nb_species) :: satol !< absolute tolerance, one value per species
    real(double_precision) :: t_stop_step

    t_stop_step = delta_t
    t = 0.d0

    do while (t.lt.t_stop_step)

      i_state = 1

      ! Adaptive absolute tolerance to avoid too high precision on very abundant
      ! species, H2 for instance. Helps running a bit faster.
      do i=1,nb_species
        satol(i) = max(a_tol, 1.d-16 * temp_abundances(i))
        if (temp_abundances(i).le.1.d-60) then
          temp_abundances(i) = 1.d-60
        end if
      enddo

      call set_work_arrays(Y=temp_abundances)

      call dlsodes(get_temporal_derivatives,nb_species,temp_abundances,t,t_stop_step,i_tol,RELATIVE_TOLERANCE,&
      satol,i_task,i_state,i_opt,rwork,lrw,iwork,liw,get_jacobian,m_f)

      ! Accumulate DLSODES step statistics (reset by set_work_arrays each call, so read now).
      tot_nst = tot_nst + int(iwork(11), 8)
      tot_nfe = tot_nfe + int(iwork(12), 8)
      tot_nje = tot_nje + int(iwork(13), 8)
      tot_nlu = tot_nlu + int(iwork(21), 8)
      tot_coldstart = tot_coldstart + 1

      ! Whenever the solver fails converging, print the reason.
      if (i_state.ne.2) then
        tot_fail = tot_fail + 1
        write(*,*) 'ISTATE = ', I_STATE
      endif

    enddo

    return
  end subroutine integrate_chemical_scheme

! =====================================================================================
! Phase III (money-plot) operator-split driver. ISOLATED from the coupled path: it is
! reached only when split_mode /= 0, and it shares nothing with the coupled integrator
! except the RHS/Jacobian (masked by the active slot range) and set_work_arrays.
! =====================================================================================

  !> Validate the split configuration and report the mask audit (per-block counts).
  subroutine split_setup_and_audit()
    implicit none
    integer :: c_lo, c_hi, d_lo, d_hi
    if (is_3_phase.ne.0) then
      write(Error_unit,'(a)') 'Error: split_mode requires the two-phase model (is_3_phase = 0).'
      call exit(41)
    endif
    if (.not.coagulation) then
      write(Error_unit,'(a)') 'Error: split_mode requires coagulation = 1 (nothing to split).'
      call exit(41)
    endif
    if (split_dt.le.0.d0) then
      write(Error_unit,'(a)') 'Error: split_mode requires split_dt > 0 [yr].'
      call exit(41)
    endif
    if (split_variant.ne.'A' .and. split_variant.ne.'B') then
      write(Error_unit,'(a)') 'Error: split_variant must be A or B.'
      call exit(41)
    endif
    if (split_order.ne.'CDC' .and. split_order.ne.'DCD') then
      write(Error_unit,'(a)') 'Error: split_order must be CDC or DCD.'
      call exit(41)
    endif
    call split_block_range(1, c_lo, c_hi)
    call split_block_range(2, d_lo, d_hi)
    ! Mask audit: the two blocks are disjoint contiguous ranges whose union is [1, nb_reactions].
    if (c_lo.ne.1 .or. d_hi.ne.nb_reactions .or. d_lo.ne.c_hi+1 .or. c_hi.lt.c_lo .or. d_hi.lt.d_lo) then
      write(Error_unit,'(a)') 'Error: split mask audit failed (blocks not a partition of the reaction set).'
      call exit(42)
    endif
    write(stdo,'(a)')       ' --- Phase III operator split ---'
    write(stdo,'(a,a,a,a)') '   variant = ', split_variant, '   order = ', split_order
    write(stdo,'(a,es12.5)') '   split_dt [yr]      = ', split_dt/YEAR
    write(stdo,'(a,i0,a,i0,a,i0)') '   chem block slots   = [', c_lo, ',', c_hi, ']  count = ', c_hi-c_lo+1
    write(stdo,'(a,i0,a,i0,a,i0)') '   dust block slots   = [', d_lo, ',', d_hi, ']  count = ', d_hi-d_lo+1
    write(stdo,'(a,i0,a,i0,a,i0,a,i0)') '   nb_chem = ', nb_chemistry_reactions, '  nb_coag_grain = ', &
          nb_coag_grain_reactions, '  nb_ice_transport = ', nb_coagulation_reactions-nb_coag_grain_reactions, &
          '  total = ', nb_reactions
  end subroutine split_setup_and_audit

  !> Slot range of a block. iblock = 1 chem, 2 dust.
  subroutine split_block_range(iblock, lo, hi)
    implicit none
    integer, intent(in)  :: iblock
    integer, intent(out) :: lo, hi
    integer :: boundary   ! last slot of the chem block
    if (split_variant.eq.'A') then
      boundary = nb_chemistry_reactions                              ! grain coag + ice -> dust
    else
      boundary = nb_chemistry_reactions + nb_coag_grain_reactions    ! ice transport only -> dust
    endif
    if (iblock.eq.1) then
      lo = 1;          hi = boundary
    else
      lo = boundary+1; hi = nb_reactions
    endif
  end subroutine split_block_range

  !> Advance one output interval by Strang macro-steps of split_dt. The state vector is
  !! inherited across sub-steps untouched (no re-initialisation); the 1e-60 floor is applied
  !! ONCE here, at the output boundary, exactly where the coupled integrator applies it.
  subroutine integrate_split(delta_t, temp_abundances, i_tol, a_tol, i_task, i_opt, m_f, i_state)
    implicit none
    real(double_precision), intent(in) :: delta_t, a_tol
    integer, intent(in) :: i_tol, i_task, i_opt, m_f
    integer, intent(out) :: i_state
    real(double_precision), dimension(nb_species), intent(inout) :: temp_abundances
    integer :: nmac, m, k
    real(double_precision) :: h
    character(len=16) :: envv
    logical, save :: split_selfcheck_done = .false.

    i_state = 2
    if (delta_t.le.0.d0) return
    nmac = nint(delta_t/split_dt)
    if (nmac.lt.1 .or. abs(dfloat(nmac)*split_dt - delta_t).gt.1.d-9*delta_t) then
      write(Error_unit,'(a,es12.5,a,es12.5,a)') 'Error: output interval ', delta_t/YEAR, &
           ' yr is not an integer multiple of split_dt = ', split_dt/YEAR, ' yr.'
      call exit(43)
    endif
    h = delta_t/dfloat(nmac)

    do k=1,nb_species          ! output-boundary floor (identical to the coupled entry floor)
      if (temp_abundances(k).le.1.d-60) temp_abundances(k) = 1.d-60
    enddo

    if (.not.split_selfcheck_done) then
      call get_environment_variable('NEMO_SPLIT_SELFCHECK', envv)
      if (len_trim(envv).gt.0) call split_selfcheck(temp_abundances)
      split_selfcheck_done = .true.
    endif

    do m=1,nmac
      if (split_order.eq.'CDC') then
        call split_substep(1, 0.5d0*h, temp_abundances, i_tol, a_tol, i_task, i_opt, m_f, i_state)
        call split_substep(2, h,       temp_abundances, i_tol, a_tol, i_task, i_opt, m_f, i_state)
        call split_substep(1, 0.5d0*h, temp_abundances, i_tol, a_tol, i_task, i_opt, m_f, i_state)
      else
        call split_substep(2, 0.5d0*h, temp_abundances, i_tol, a_tol, i_task, i_opt, m_f, i_state)
        call split_substep(1, h,       temp_abundances, i_tol, a_tol, i_task, i_opt, m_f, i_state)
        call split_substep(2, 0.5d0*h, temp_abundances, i_tol, a_tol, i_task, i_opt, m_f, i_state)
      endif
      tot_macro = tot_macro + 1
    enddo

    active_lo = 1              ! restore the full range (coupled state) on exit
    active_hi = nb_reactions
  end subroutine integrate_split

  !> One sub-step: a FRESH DLSODES solve (ISTATE=1, work arrays rebuilt) of the masked
  !! system over [0, h]. Solver history is never carried across a mask switch; the
  !! abundances are. No floor here (sec. 2c). Restarts after a failure also cold-start.
  subroutine split_substep(iblock, h, Y, i_tol, a_tol, i_task, i_opt, m_f, i_state)
    implicit none
    integer, intent(in) :: iblock, i_tol, i_task, i_opt, m_f
    real(double_precision), intent(in) :: h, a_tol
    real(double_precision), dimension(nb_species), intent(inout) :: Y
    integer, intent(out) :: i_state
    real(double_precision) :: t
    real(double_precision), dimension(nb_species) :: satol
    integer :: k

    call split_block_range(iblock, active_lo, active_hi)
    t = 0.d0
    do while (t.lt.h)
      i_state = 1
      do k=1,nb_species
        satol(k) = max(a_tol, 1.d-16 * Y(k))
      enddo
      call set_work_arrays(Y=Y)
      call dlsodes(get_temporal_derivatives,nb_species,Y,t,h,i_tol,RELATIVE_TOLERANCE,&
      satol,i_task,i_state,i_opt,rwork,lrw,iwork,liw,get_jacobian,m_f)
      tot_nst = tot_nst + int(iwork(11), 8)
      tot_nfe = tot_nfe + int(iwork(12), 8)
      tot_nje = tot_nje + int(iwork(13), 8)
      tot_nlu = tot_nlu + int(iwork(21), 8)
      tot_coldstart = tot_coldstart + 1
      if (i_state.ne.2) then
        tot_fail = tot_fail + 1
        write(*,*) 'ISTATE = ', i_state, ' (split block ', iblock, ')'
      endif
    enddo
  end subroutine split_substep

  !> Debug self-check (NEMO_SPLIT_SELFCHECK=1), split mode only, once per run:
  !!  (1) additivity  f_chem + f_dust == f_full  and  J_chem + J_dust == J_full (every column);
  !!  (2) dust-block analytic Jacobian vs central FD of the masked dust RHS. The dust block is
  !!      bilinear, so the central difference is exact up to round-off, and with chemistry masked
  !!      the tiny coagulation entries are resolvable (they are not against the full RHS).
  subroutine split_selfcheck(Y)
    implicit none
    real(double_precision), dimension(nb_species), intent(in) :: Y
    real(double_precision), dimension(nb_species) :: fF, fC, fD, jF, jC, jD, fp, fm, Yp, jfd
    real(double_precision) :: dum_ian(3), dum_jan(3), h, e_add_f, e_add_j, e_fd, scale, colmax
    integer :: j, i, lo, hi, n_extra, n_missing, n_sig, n_unres
    real(double_precision) :: fscale
    integer :: save_lo, save_hi
    save_lo = active_lo; save_hi = active_hi
    call set_constant_rates()
    active_lo = 1; active_hi = nb_reactions;   call get_temporal_derivatives(nb_species, 0.d0, Y, fF)
    call split_block_range(1, lo, hi); active_lo = lo; active_hi = hi; call get_temporal_derivatives(nb_species, 0.d0, Y, fC)
    call split_block_range(2, lo, hi); active_lo = lo; active_hi = hi; call get_temporal_derivatives(nb_species, 0.d0, Y, fD)
    scale = maxval(abs(fF))
    e_add_f = maxval(abs(fC + fD - fF)) / scale
    e_add_j = 0.d0; e_fd = 0.d0; n_extra = 0; n_missing = 0; n_sig = 0; n_unres = 0
    do j=1,nb_species
      active_lo = 1; active_hi = nb_reactions
      call get_temporal_derivatives(nb_species, 0.d0, Y, fp)          ! refresh dependant rates at Y
      call get_jacobian(3, 0.d0, Y, j, dum_ian, dum_jan, jF)
      call split_block_range(1, lo, hi); active_lo = lo; active_hi = hi
      call get_temporal_derivatives(nb_species, 0.d0, Y, fp)
      call get_jacobian(3, 0.d0, Y, j, dum_ian, dum_jan, jC)
      call split_block_range(2, lo, hi); active_lo = lo; active_hi = hi
      call get_jacobian(3, 0.d0, Y, j, dum_ian, dum_jan, jD)
      colmax = max(maxval(abs(jF)), tiny(1.d0))
      e_add_j = max(e_add_j, maxval(abs(jC + jD - jF)) / colmax)
      ! central FD of the masked (dust) RHS in column j
      h = 1.d-4 * abs(Y(j)); if (h.eq.0.d0) cycle
      Yp = Y; Yp(j) = Y(j) + h; call get_temporal_derivatives(nb_species, 0.d0, Yp, fp)
      Yp = Y; Yp(j) = Y(j) - h; call get_temporal_derivatives(nb_species, 0.d0, Yp, fm)
      jfd = (fp - fm) / (2.d0*h)
      ! An entry is RESOLVABLE by the FD only if its induced change 2h|J| clears the round-off
      ! of f_i by a wide margin (1e-8 relative); below that the FD is 0 or noise by construction.
      do i=1,nb_species
        fscale = max(abs(fp(i)), abs(fm(i)), tiny(1.d0))
        if (2.d0*h*max(abs(jD(i)),abs(jfd(i))).gt.1.d-8*fscale) then
          n_sig = n_sig + 1
          e_fd = max(e_fd, abs(jD(i)-jfd(i)) / max(abs(jD(i)),abs(jfd(i))))
          if (jD(i).eq.0.d0) n_missing = n_missing + 1
          if (jfd(i).eq.0.d0) n_extra = n_extra + 1
        elseif (jD(i).ne.0.d0) then
          n_unres = n_unres + 1
        endif
      enddo
    enddo
    active_lo = save_lo; active_hi = save_hi
    write(stdo,'(a)') ' --- split self-check (NEMO_SPLIT_SELFCHECK) ---'
    write(stdo,'(a,es10.3)') '   RHS additivity  max|f_C+f_D-f_full|/max|f_full|     = ', e_add_f
    write(stdo,'(a,es10.3)') '   Jac additivity  max_col max|J_C+J_D-J_full|/max|J_full| = ', e_add_j
    write(stdo,'(a,i0,a,es10.3)') '   dust Jac vs FD  (', n_sig, ' significant entries) max rel diff = ', e_fd
    write(stdo,'(a,i0,a,i0)') '   dust Jac entries missing (FD>0, J=0) = ', n_missing, '   extra (J>0, FD=0) = ', n_extra
    write(stdo,'(a,i0)') '   dust Jac entries below FD resolution (not tested)   = ', n_unres
  end subroutine split_selfcheck

END program nmgc


!-------------------------------------------------------------------------
!              GET ARGUMENT EITHER FROM ARGV OR FROM STDI
!-------------------------------------------------------------------------
subroutine ggetarg(iarg,buffer,fromstdi)
  use global_variables
  implicit none
  character*100 :: buffer
  integer :: iarg
  logical :: fromstdi
  if (fromstdi) then
     read(ffli,*,end=103) buffer
  else
     call getarg(iarg,buffer)
  endif
  return
103 continue
  buffer='enter'
  return
end subroutine ggetarg


!-------------------------------------------------------------------------
!                    INTERPRET COMMAND-LINE OPTIONS
!-------------------------------------------------------------------------
subroutine interpet_command_line(gotit,fromstdi,quit)
  use global_variables
  implicit none
  character*100 :: buffer
  integer :: iarg,numarg
  logical :: gotit,quit,fromstdi
  gotit = .false.
  quit  = .false.

  if (fromstdi) then
    numarg = 1000
  else
    numarg = iargc()
  endif
  iarg = 1
  do while(iarg.le.numarg)
    call ggetarg(iarg,buffer,fromstdi)
    iarg = iarg+1

    if (buffer(1:3).eq.'run') then
      do_chemical_scheme = .true.
      gotit = .true.
    elseif (buffer(1:7).eq.'outputs') then
      do_outputs = .true.
      gotit = .true.
    else
      write(stdo,*) 'ERROR: Could not recognize command line option ',trim(buffer)
      stop
    end if
  end do

  if (.not.gotit) then
    call write_banner()
    write(stdo,*) 'Please, use one of these actions:'
    write(stdo,*) '  run        : Integrate the evolution of the chemical scheme'
    write(stdo,*) '  outputs    : Read binary outputs to convert into ASCII format'
    quit = .true.
  endif

end subroutine interpet_command_line


!-------------------------------------------------------------------------
!                  WRITE THE BANNER ON THE SCREEN
!-------------------------------------------------------------------------
subroutine write_banner()
  use global_variables
  implicit none
  write(stdo,*) ' '
  write(stdo,*) '================================================================'
  write(stdo,*) '        0D MULTI-GRAIN GAS-GRAIN CODE  --  skeleton stage        '
  write(stdo,*) '            extracted from NMGC-2.0 (Gavino et al.)              '
  write(stdo,*) '================================================================'
  write(stdo,*) ' '
  call flush(stdo)
end subroutine write_banner
