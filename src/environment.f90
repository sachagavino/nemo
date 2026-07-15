!******************************************************************************
! MODULE: environment
!******************************************************************************
!
! DESCRIPTION:
!> @brief Physical environment of the 0D box: gas properties (density, Av,
!! temperature) and the size -> property maps of the dust population.
!!
!! This module replaces the former 'structure' module of nmgc-2.0, from which
!! every 1D construct (spatial grid, species diffusion, Crank-Nicholson,
!! 1D_static.dat / structure_evolution.dat tables) has been removed.
!!
!! HOOKS (deliberately left in place, nothing implemented):
!!   * get_structure_properties is still reached through a procedure pointer
!!     (see global_variables), so a time-dependent 0D structure (collapsing
!!     core, parcel trajectory) can be re-added as one extra case without
!!     touching any caller.
!!   * the grain temperature prescriptions below are the place where Td(a) and
!!     T_CR,peak(a) will become explicit, named, switchable prescriptions
!!     evaluated once on the size grid at init (stage 3).
!
!******************************************************************************

module environment

use iso_fortran_env
use numerical_types
use global_variables

implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Properties of the 0D box (Av, gas temperature, gas density).
!! Here, everything is constant in time.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine get_structure_properties_fixed(time, Av, density, gas_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: time !<[in] Current time of the simulation [s]

  ! Outputs
  real(double_precision), intent(out) :: Av !<[out] Visual extinction [mag]
  real(double_precision), intent(out) :: gas_temperature !<[out] gas temperature [K]
  real(double_precision), intent(out) :: density !<[out] gas density [part/cm^3]

  !----------------------------------------------------------------------------

  density = initial_gas_density
  av = initial_visual_extinction
  gas_temperature = initial_gas_temperature

  return
end subroutine get_structure_properties_fixed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Grain temperature, assumed equal to the gas temperature.
!! Same value for every size bin.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine get_grain_temperature_gas(time, gas_temperature, av, grain_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
  real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
  real(double_precision), intent(in) :: av !<[in] visual extinction [mag]

  ! Outputs
  real(double_precision), dimension(:), intent(out) :: grain_temperature !<[out] dim(nb_grains) [K]
  !----------------------------------------------------------------------------

  grain_temperature(1:nb_grains) = gas_temperature

  return
end subroutine get_grain_temperature_gas

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Grain temperature fixed to initial_dust_temperature, identical for
!! every size bin. Only meaningful in single-grain mode.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine get_grain_temperature_fixed(time, gas_temperature, av, grain_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
  real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
  real(double_precision), intent(in) :: av !<[in] visual extinction [mag]

  ! Outputs
  real(double_precision), dimension(:), intent(out) :: grain_temperature !<[out] dim(nb_grains) [K]
  !----------------------------------------------------------------------------

  if (multi_grain.eq.1) then
    write(error_unit,*) 'Please check parameters.in:'
    write(error_unit,*) 'if multi_grain = 1, grain_temperature_type = fixed is not valid.'
    call exit(1)
  else
    grain_temperature(1:nb_grains) = initial_dust_temperature
  endif

  return
end subroutine get_grain_temperature_fixed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Grain temperature set per size bin.
!!
!! In this extraction stage the per-bin values are still whatever was loaded
!! into grain_temp(:) at init (i.e. the 3rd column of 0D_grain_sizes.in), so
!! that the skeleton reproduces nmgc-2.0 exactly.
!!
!! STAGE 3: this becomes a real physical statement, Td = Td(a), evaluated once
!! on the size grid from a named, switchable prescription (e.g. the ad-hoc
!! Td ~ a^(-1/6) scaling currently hard-coded in the time loop of main.f90 of
!! nmgc-2.0, promoted here and given an explicit reference radius).
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine get_grain_temperature_fixed_to_dust_size(time, gas_temperature, av, grain_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
  real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
  real(double_precision), intent(in) :: av !<[in] visual extinction [mag]

  ! Outputs
  real(double_precision), dimension(:), intent(out) :: grain_temperature !<[out] dim(nb_grains) [K]

  ! Locals
  integer :: i
  !----------------------------------------------------------------------------

  if (multi_grain.eq.0) then
    write(error_unit,*) 'Please check parameters.in:'
    write(error_unit,*) 'if multi_grain = 0, grain_temperature_type = fixed_to_dust_size is not valid.'
    call exit(1)
  endif

  do i=1,nb_grains
    grain_temperature(i) = grain_temp(i)
  enddo

  return
end subroutine get_grain_temperature_fixed_to_dust_size

end module environment
