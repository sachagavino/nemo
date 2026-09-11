program r3_charge_check
  use global_variables
  use gasgrain
  implicit none
  integer :: k, r
  character(len=11) :: g0,gm
  real(double_precision) :: rf, rd
  call init_gasgrain()
  write(*,'(a)') 'CHARGE-RATE a^2 CHECK'
  write(*,'(a)') ' bin  a[cm]      k_form/a^2       k_dest/a^2'
  do k=1,nb_grains
    write(g0,'(a,i2.2)') 'GRAIN',k
    write(gm,'(a,i2.2,a)') 'GRAIN',k,'-'
    rf=0; rd=0
    do r=1,nb_reactions
      if (trim(REACTION_COMPOUNDS_NAMES(1,r))==trim(g0).and.trim(REACTION_COMPOUNDS_NAMES(2,r))=='e-') rf=reaction_rates(r)
      if (rd==0.d0.and.trim(REACTION_COMPOUNDS_NAMES(1,r))==trim(gm)) rd=reaction_rates(r)
    enddo
    write(*,'(i4,3es16.6)') k, grain_radii(k), rf/grain_radii(k)**2, rd/grain_radii(k)**2
  enddo
end program
