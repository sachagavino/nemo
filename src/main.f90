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
  use utilities
  use gasgrain
  use outputs

  implicit none

  ! Parameters for DLSODES. RTOL is the RELATIVE_TOLERANCE parameter in global_variables.f90
  integer :: itol = 2 !< ITOL = 1 or 2 according as ATOL is a scalar or array.
  integer :: itask = 1 !< ITASK = 1 for normal computation of output values of Y at t = TOUT.
  integer :: istate = 1 !< ISTATE = integer flag (input and output). Set ISTATE = 1.
  integer :: iopt = 1 !< IOPT = 1 to indicate optional inputs are used.
  integer :: mf = 121 !< method flag. 121: stiff (BDF) method, user-supplied sparse Jacobian
  real(double_precision) :: atol = 1.d-99 !< absolute tolerance parameter

  real(double_precision) :: output_timestep !< Timestep to reach the next output time [s]

  real(double_precision), dimension(:), allocatable :: output_times !< [s] dim(NB_OUTPUTS)
  real(double_precision) :: output_step !< [s] linear or log spacing of the output times
  integer :: output_idx !< index of the output time loop

  real(double_precision) :: code_start_time, code_current_time, code_elapsed_time
  real(double_precision) :: remaining_time !< estimated remaining time [s]

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

      ! Column densities of the self-shielded species, from Av
      NH    = actual_av/AV_NH_ratio * abundances(indH)
      NH2   = actual_av/AV_NH_ratio * abundances(indH2)
      NN2   = actual_av/AV_NH_ratio * abundances(indN2)
      NCO   = actual_av/AV_NH_ratio * abundances(indCO)
      NH2O  = actual_av/AV_NH_ratio * abundances(indH2O)
      NCH   = actual_av/AV_NH_ratio * abundances(indCH)
      NCH3  = actual_av/AV_NH_ratio * abundances(indCH3)
      NH2CO = actual_av/AV_NH_ratio * abundances(indH2CO)
      NCO2  = actual_av/AV_NH_ratio * abundances(indCO2)
      NN2O  = actual_av/AV_NH_ratio * abundances(indN2O)
      NCH4  = actual_av/AV_NH_ratio * abundances(indCH4)
      NOH   = actual_av/AV_NH_ratio * abundances(indOH)
      NHCO  = actual_av/AV_NH_ratio * abundances(indHCO)
      NCN   = actual_av/AV_NH_ratio * abundances(indCN)
      NHCN  = actual_av/AV_NH_ratio * abundances(indHCN)
      NHNC  = actual_av/AV_NH_ratio * abundances(indHNC)
      NNH   = actual_av/AV_NH_ratio * abundances(indNH)
      NNH2  = actual_av/AV_NH_ratio * abundances(indNH2)
      NNH3  = actual_av/AV_NH_ratio * abundances(indNH3)

      call integrate_chemical_scheme(delta_t=output_timestep, temp_abundances=abundances(1:nb_species), & ! Inputs
      i_tol=itol, a_tol=atol, i_task=itask, i_opt=iopt, m_f=mf, & ! Inputs
      i_state=istate) ! Output

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
      call write_current_dust(index=output_idx)

      first_step_done = .true.
    enddo

    call write_abundances('abundances.tmp')

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

      ! Whenever the solver fails converging, print the reason.
      if (i_state.ne.2) then
        write(*,*) 'ISTATE = ', I_STATE
      endif

    enddo

    return
  end subroutine integrate_chemical_scheme

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
