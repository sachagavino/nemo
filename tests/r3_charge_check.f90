program r3_charge_check
  use global_variables
  use gasgrain
  implicit none
  integer :: k, r
  character(len=11) :: gname
  call init_gasgrain()
  write(*,'(a)') 'per-bin electron-sticking rate  GRAIN0k + e- -> GRAINk-'
  write(*,'(a)') ' bin   a[cm]         reaction_rate        rate/a^2'
  do k=1,nb_grains
    write(gname,'(a,i2.2)') 'GRAIN', k
    do r=1,nb_reactions
      if (trim(REACTION_COMPOUNDS_NAMES(1,r))==trim(gname) .and. trim(REACTION_COMPOUNDS_NAMES(2,r))=='e-') then
        write(*,'(i4,3es18.6)') k, grain_radii(k), reaction_rates(r), reaction_rates(r)/grain_radii(k)**2
      endif
    enddo
  enddo
end program r3_charge_check
