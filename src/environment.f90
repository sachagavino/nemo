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
!> @brief Grain temperature set per size bin.
!!
!! In this extraction stage the per-bin values are still whatever was loaded
!! into grain_temp(:) at init (tabulated path: 3rd column of dust_grid_table.in), so
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

  do i=1,nb_grains
    grain_temperature(i) = grain_temp(i)
  enddo

  return
end subroutine get_grain_temperature_fixed_to_dust_size

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @brief Named size->temperature prescriptions, evaluated ONCE on the grid at
!! init (constant per bin, fixed in time for now). This is what makes
!! fixed_to_dust_size on a DERIVED grid a real physical statement Td = Td(a)
!! rather than a table read. On a tabulated grid the per-bin temperatures come
!! from dust_grid_table.in instead and this evaluator is not called.
!!
!! Default prescription: the ad-hoc Td ~ a^(-1/6) scaling that nmgc-2.0 buried
!! in the main.f90 time loop, promoted here and anchored to an explicit
!! reference: Td(a) = reference_dust_temperature * (a/reference_grain_radius)^(-1/6).
!! T_CR,peak(a) is uniform for now (= cr_peak_grain_temp); the hook is here to
!! make it size-dependent later.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function td_of_a(a) result(td)
  implicit none
  real(double_precision), intent(in) :: a !<[in] grain radius [cm]
  real(double_precision) :: td !< dust temperature [K]
  td = reference_dust_temperature * (a/reference_grain_radius)**(-1.d0/6.d0)
end function td_of_a

pure function tcr_peak_of_a(a) result(tcr)
  implicit none
  real(double_precision), intent(in) :: a !<[in] grain radius [cm]
  real(double_precision) :: tcr !< CR peak grain temperature [K]
  tcr = cr_peak_grain_temp   ! uniform for now; size-dependent hook
end function tcr_peak_of_a

subroutine evaluate_dust_temperature_prescriptions()
  implicit none
  integer :: i
  do i=1,nb_grains
    grain_temp(i)              = td_of_a(grain_radii(i))
    CR_PEAK_GRAIN_TEMP_all(i)  = tcr_peak_of_a(grain_radii(i))
  enddo
  return
end subroutine evaluate_dust_temperature_prescriptions

end module environment
