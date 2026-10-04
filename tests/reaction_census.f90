! ===========================================================================
! tests/reaction_census.f90 -- per-type reaction and species census.
!
! Calls init_gasgrain() on the run-directory config and prints, straight from the
! reaction generator's own counters and the state vector, the counts quoted in the
! paper's dimensionality paragraph:
!   - species: gas / ice / grain (and total)
!   - reactions: gas-chemistry / surface-chemistry (split of nb_chemistry_reactions)
!                coagulation (grain) / ice-transport  (the generated COAGULATION_TYPE
!                block), and the totals.
! No integration. Run once per N (set by a_max in parameters.in); quote from the run.
! ===========================================================================
program reaction_census
use iso_fortran_env, only: output_unit
use global_variables
use gasgrain
implicit none

integer :: i, r, s, idx
integer :: ngas, nice, ngrain
integer :: nc_grain, nc_ice, nchem_gas, nchem_surf
logical :: is_coag, touches_surface
character(len=len(species_name(1))) :: nm

call init_gasgrain()

! ---- species census (by name) ----
ngas = 0 ; nice = 0 ; ngrain = 0
do i = 1, nb_species
  nm = species_name(i)
  if (is_grain_name(nm)) then
    ngrain = ngrain + 1
  else if (is_ice_name(nm)) then
    nice = nice + 1
  else
    ngas = ngas + 1
  endif
enddo

! ---- reaction census (by type / species content) ----
nc_grain = 0 ; nc_ice = 0 ; nchem_gas = 0 ; nchem_surf = 0
do r = 1, nb_reactions
  is_coag = (REACTION_TYPE(r) == COAGULATION_TYPE)
  if (is_coag) then
    ! grain-only vs ice-transport by reactant-1 (GRAIN* vs J*/K*)
    idx = REACTION_COMPOUNDS_ID(1, r)
    if (valid(idx)) then
      if (is_grain_name(species_name(idx))) then
        nc_grain = nc_grain + 1
      else if (is_ice_name(species_name(idx))) then
        nc_ice = nc_ice + 1
      endif
    endif
  else
    ! chemistry: surface if any compound is a J/ice or GRAIN species, else gas
    touches_surface = .false.
    do s = 1, MAX_COMPOUNDS
      idx = REACTION_COMPOUNDS_ID(s, r)
      if (.not. valid(idx)) cycle
      if (is_ice_name(species_name(idx)) .or. is_grain_name(species_name(idx))) then
        touches_surface = .true. ; exit
      endif
    enddo
    if (touches_surface) then
      nchem_surf = nchem_surf + 1
    else
      nchem_gas = nchem_gas + 1
    endif
  endif
enddo

write(output_unit,'(a)') '================ reaction / species census ================'
write(output_unit,'(a,i0)')  '  nb_grains (N)                 = ', nb_grains
write(output_unit,'(a)') '  --- species ---'
write(output_unit,'(a,i0)')  '    gas                         = ', ngas
write(output_unit,'(a,i0)')  '    ice (per bin x N)           = ', nice
write(output_unit,'(a,i0)')  '    grain (charge states x N)   = ', ngrain
write(output_unit,'(a,i0)')  '    TOTAL nb_species            = ', nb_species
write(output_unit,'(a)') '  --- reactions ---'
write(output_unit,'(a,i0)')  '    gas chemistry               = ', nchem_gas
write(output_unit,'(a,i0)')  '    surface chemistry           = ', nchem_surf
write(output_unit,'(a,i0)')  '    chemistry subtotal (gas+surf)= ', nchem_gas + nchem_surf
write(output_unit,'(a,i0)')  '      [nb_chemistry_reactions]   = ', nb_chemistry_reactions
write(output_unit,'(a,i0)')  '    coagulation (grain pairs)    = ', nc_grain
write(output_unit,'(a,i0)')  '      [nb_coag_grain_reactions]  = ', nb_coag_grain_reactions
write(output_unit,'(a,i0)')  '    ice transport                = ', nc_ice
write(output_unit,'(a,i0)')  '      [nb_ice_transport_reactions]= ', nb_ice_transport_reactions
write(output_unit,'(a,i0)')  '    coag+ice block subtotal      = ', nc_grain + nc_ice
write(output_unit,'(a,i0)')  '    TOTAL nb_reactions           = ', nb_reactions
write(output_unit,'(a)') '==========================================================='

contains

  logical function valid(k)
    integer, intent(in) :: k
    valid = (k >= 1 .and. k <= nb_species)
  end function valid

  logical function is_grain_name(nm_in)
    character(len=*), intent(in) :: nm_in
    is_grain_name = .false.
    if (len_trim(nm_in) >= 5) then
      if (nm_in(1:5) == 'GRAIN') is_grain_name = .true.
    endif
  end function is_grain_name

  logical function is_ice_name(nm_in)
    character(len=*), intent(in) :: nm_in
    is_ice_name = (nm_in(1:1) == 'J' .or. nm_in(1:1) == 'K')
  end function is_ice_name

end program reaction_census
