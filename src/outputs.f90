!******************************************************************************
! MODULE: outputs
!******************************************************************************
!
! DESCRIPTION:
!> @brief Convert the binary output of a 0D run (abundances.out) into ASCII
!! files: one file per species in ab/, one file per surface/mantle species in
!! ml/ (monolayers), and the physical structure history in struct/.
!!
!! HOOK (coagulation): the dust diagnostics required from day one -- dust mass
!! sum_k m_k Y_k, grain number budget, and total ice per species summed over
!! bins -- will be written from here, once the dust bins are part of the state
!! vector.
!
!******************************************************************************

module outputs

use iso_fortran_env
use numerical_types
use utilities
use gasgrain

implicit none

contains

subroutine get_outputs()
  implicit none
  ! Locals
  character(len=80) :: filename_output
  integer :: species, output
  logical :: isDefined

  character(len=80) :: output_format, output_format1, output_format2
  character(len=512) :: sys_cmd

  real(double_precision), dimension(:,:), allocatable :: abundances_out !< dim(nb_outputs, nb_species)
  real(double_precision), dimension(:), allocatable :: time !< [s]
  real(double_precision), dimension(:), allocatable :: gas_temperature_out !< [K]
  real(double_precision), dimension(:,:), allocatable :: dust_temperature_out !< dim(nb_grains, nb_outputs) [K]
  real(double_precision), dimension(:), allocatable :: density !< [part/cm^3]
  real(double_precision), dimension(:), allocatable :: visual_extinction_out !< [mag]
  real(double_precision), dimension(:), allocatable :: x_rate !< X ionisation rate [s-1]
  character(2) :: c_i
  integer :: i, ic_i
  integer :: ios !< IOSTAT flag for detecting early end-of-file

  allocate(time(nb_outputs))
  allocate(gas_temperature_out(nb_outputs))
  allocate(dust_temperature_out(nb_grains, nb_outputs))
  allocate(density(nb_outputs))
  allocate(visual_extinction_out(nb_outputs))
  allocate(x_rate(nb_outputs))
  allocate(abundances_out(nb_outputs, nb_species))

  write(*,'(a)', advance='no') 'Reading unformatted outputs...'
  ! All timesteps are stored sequentially. The file may contain fewer timesteps
  ! than nb_outputs if the run was interrupted.
  open(10, file='abundances.out', status='old', form='unformatted')
  do output=1,nb_outputs
    read(10, iostat=ios) time(output)
    if (ios /= 0) then
      nb_outputs = output - 1
      exit
    endif
    read(10) gas_temperature_out(output), dust_temperature_out(1:nb_grains, output), &
             density(output), visual_extinction_out(output), x_rate(output)
    read(10) abundances_out(output, 1:nb_species)
  enddo
  close(10)
  write(*,'(a,a)') achar(13), 'Reading unformatted outputs... Done'
  write(*,'(a,i0,a)') 'Note: ', nb_outputs, ' timestep(s) found in abundances.out'

  !####################################################
  ! One file per species, one line per output time
  !####################################################
  inquire(file='ab', exist=isDefined)
  if (.not.isDefined) call system("mkdir ab")

  inquire(file='ml', exist=isDefined)
  if (.not.isDefined) call system("mkdir ml")

  write(filename_output, '(a)') 'ab/gas_phase.ab'
  open(11, file=filename_output)
  write(11,'(a)')'! species_name  abundance@last_output. &
                & Note:- outputs are in descending order w.r.t. second column'
  write(output_format1, *) '(A,4X,es13.6e2)'
  write(output_format, *) '(es10.3e2,es13.6e2)'

  write(*,'(a)', advance='no') 'Writing *.ab ASCII files in ab/...'
  do species=1, nb_species
    write(filename_output, '(a,a,a)') 'ab/', trim(species_name(species)), '.ab'
    open(10, file=filename_output)
    write(10,'(a)') '! time [year] ; abundance (relative to H) [number ratio]'

    do output=1, nb_outputs
      write(10,output_format) time(output)/YEAR, abundances_out(output, species)
    enddo
    close(10)

    if (species_name(species)(1:1) /= 'J' .and. species_name(species)(1:1) /= 'K') &
      write(11,output_format1) species_name(species), abundances_out(nb_outputs, species)
  enddo
  close(11)

  ! Sort ab/gas_phase.ab from high to low abundance
  write(filename_output, '(a)') 'ab/gas_phase.ab'
  sys_cmd = 'head -1 '//trim(filename_output)//' > ab/ab.tmp'
  call system(sys_cmd)
  sys_cmd = 'tail -n+2 '//trim(filename_output)//' | sort  -k2,2 -gr >>  ab/ab.tmp'
  call system(sys_cmd)
  sys_cmd = 'mv ab/ab.tmp '//trim(filename_output)
  call system(sys_cmd)
  write(*,'(a,a)') achar(13), 'Writing output files in ab/... Done'

  !####################################################
  ! Monolayers, one file per surface/mantle species
  !####################################################
  write(*,'(a)', advance='yes') 'Writing mono layer data for grain species in *.ml ASCII files in the folder ml/...'

  do i=1,nb_grains
    write(filename_output, '(a,i2.2,a)') 'ml/grain',i,'.ml'
    open(100+i, file=filename_output)
    write(100+i,'(a)')'! species_name  mono-layers@last_output  abundance@last_output. &
                      & Note:- outputs are in descending order w.r.t. mono-layers, i.e. second column.'
  enddo

  write(output_format, *) '(es10.3e2,es13.6e2)'
  write(output_format2, *)'(A,2X,es13.3e2,13X,es13.6e2)'

  do species=1, nb_species
    if (species_name(species)(1:1) == 'J' .or. species_name(species)(1:1) == 'K') then
      write(filename_output, '(a,a,a)') 'ml/', trim(species_name(species)), '.ml'
      open(10, file=filename_output)
      write(10,'(a)') '! time [year] ; number of mono layers'
      c_i = species_name(species)(2:3)
      read(c_i,'(I2)') ic_i
      do output=1, nb_outputs
        write(10,output_format) time(output)/YEAR, &
        abundances_out(output, species) * GTODN(ic_i) / nb_sites_per_grain(ic_i)   ! abundance -> mono layers
      enddo
      close(10)

      write(100+ic_i,output_format2) species_name(species), &
      abundances_out(nb_outputs, species) * GTODN(ic_i) / nb_sites_per_grain(ic_i), &
      abundances_out(nb_outputs, species)
    endif
  enddo

  write(*,'(a,a)') achar(13), 'Writing output files in the folder ml/... Done'

  do i=1,nb_grains
    write(filename_output, '(a,i2.2,a)') 'ml/grain',i,'.ml'
    sys_cmd = 'head -1 '//trim(filename_output)//' > ml/ml.tmp'
    call system(sys_cmd)
    sys_cmd = 'tail -n+2 '//trim(filename_output)//' | sort  -k2,2 -gr >>  ml/ml.tmp'
    call system(sys_cmd)
    sys_cmd = 'mv ml/ml.tmp '//trim(filename_output)
    call system(sys_cmd)
    ! sum the monolayers of all species on that grain -> ml/sum_all_sp.ml
    sys_cmd = "awk 'END { print "//'"'//trim(filename_output(4:10))&
    //'"'//" ,s}{s+=$2}' "//trim(filename_output)//'>> ml/sum_all_sp.ml'
    call system(sys_cmd)
    close(100+i)
  enddo

  !####################################################
  ! Physical structure history
  !####################################################
  inquire(file='struct', exist=isDefined)
  if (.not.isDefined) call system("mkdir struct")
  call system("rm -f struct/*.struct")

  write(*,'(a)', advance='no') 'Writing *.struct ASCII files in struct/...'
  write(filename_output, '(a)') 'struct/output.struct'
  open(10, file=filename_output)
  write(10,'(a)') '! time     ; gas temperature ; dust temperature (bin 1)&
                  & ; H density  ; visual extinction ; x ionization rate'
  write(10,'(a)') '!  [year]  ;       [K]       ;         [K]     &
                  & ; [part/cm^3] ;           [mag]   ;       [s-1] '

  do output=1, nb_outputs
    write(10,'(6(es10.3e2," "))') time(output)/YEAR, gas_temperature_out(output), dust_temperature_out(1, output), &
          density(output), visual_extinction_out(output), x_rate(output)
  enddo
  close(10)
  write(*,'(a,a)') achar(13), 'Writing structure output files in struct/... Done'

end subroutine get_outputs

end module outputs
