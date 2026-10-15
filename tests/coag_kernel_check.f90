! ===========================================================================
! tests/coag_kernel_check.f90
!
! Absolute check of the free-molecular Brownian coagulation kernel, to pin the
! prefactor convention (the mean-vs-RMS relative-speed ambiguity). Initialises the
! dust-only fixture with coagulation_kernel=brownian and compares coag_kernel(i,j)
! for two specific pairs against values hand-computed OFFLINE for this exact grid
! and temperature -- so a wrong prefactor (e.g. RMS speed, or a dropped/extra pi or
! factor of 2) fails here even though it would still "look Brownian".
!
! Fixture (tests/coag_kernel_check.sh sets it): a1 = 1e-6 cm, mass_ratio = 2,
! grain_density = 3 g/cm^3, gas T = 27.41 K. Then, with mu_ij = m_i m_j/(m_i+m_j),
!   K_ij   = (a_i+a_j)^2 * sqrt( 8 pi k_B T / mu_ij )        [free molecular, MEAN speed]
!   K(1,1) = 1/2 * K_11^phys  (self-pair counted once)
! Hand-computed (k_B = 1.3806488e-16 cgs):
!   K(1,1) = 2.4606855503e-10 cm^3/s
!   K(1,2) = 5.4418089399e-10 cm^3/s
! (The RMS-speed form would give K(1,2) ~ 5.907e-10, 8.5% high -- caught at 0.5%.)
!
! Build/run via tests/coag_kernel_check.sh. Expected: RESULT: PASS.
! ===========================================================================
program coag_kernel_check
  use global_variables
  use gasgrain
  use dust_evolution
  implicit none
  real(double_precision) :: k11, k12, e11, e12, r11, r12, tol
  integer :: nbug

  call init_gasgrain()          ! reads the Brownian fixture from the working directory

  e11 = 2.4606855503d-10        ! hand-computed absolute references (see header)
  e12 = 5.4418089399d-10
  tol = 5.d-3                    ! 0.5% -- tight enough to reject the RMS-speed prefactor

  if (.not. coagulation) then
    write(*,'(a)') ' RESULT: FAIL (coagulation is off -- fixture misconfigured)'
    call exit(1)
  endif
  if (trim(coagulation_kernel) /= 'brownian') then
    write(*,'(3a)') ' RESULT: FAIL (coagulation_kernel = "', trim(coagulation_kernel), '", expected brownian)'
    call exit(1)
  endif

  k11 = coag_kernel(1,1)
  k12 = coag_kernel(1,2)
  r11 = abs(k11 - e11)/e11
  r12 = abs(k12 - e12)/e12
  nbug = 0
  if (r11 > tol) nbug = nbug + 1
  if (r12 > tol) nbug = nbug + 1

  write(*,'(a)') '=================================================================='
  write(*,'(a)') ' Brownian coagulation kernel absolute check (free molecular, mean speed)'
  write(*,'(a)') '=================================================================='
  write(*,'(a,f8.3,a)')                      ' gas temperature      = ', gas_temperature, ' K'
  write(*,'(a,es16.8,a,es16.8,a,es9.2)') ' K(1,1) self : got ', k11, '  expect ', e11, '  relerr ', r11
  write(*,'(a,es16.8,a,es16.8,a,es9.2)') ' K(1,2) cross: got ', k12, '  expect ', e12, '  relerr ', r12
  write(*,'(a,es9.2)')                       ' tolerance            = ', tol
  if (nbug == 0) then
    write(*,'(a)') ' RESULT: PASS'
  else
    write(*,'(a)') ' RESULT: FAIL'
    call exit(1)
  endif
end program coag_kernel_check
