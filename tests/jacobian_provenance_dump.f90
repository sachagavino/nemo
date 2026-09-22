! ===========================================================================
! tests/jacobian_provenance_dump.f90 -- Phase III (docs/PhaseIII_jacobian_figure_brief.md)
!
! READ-ONLY dump of the symbolic Jacobian sparsity pattern the solver already
! builds (IA_SYM/JA_SYM, via build_symbolic_sparsity in init_gasgrain). It does
! NOT modify any production path and does NOT integrate: init only.
!
! For every nonzero pattern cell (row,col) it computes a bitmask of ALL
! contributing mechanisms (overlaps are the point) and writes:
!   jac_pattern.tsv  : row  col  mech_mask     (one line per nonzero)
!   jac_species.tsv  : index  name  class  bin  charge
!
! Classifier (species-content, per the design-thread disposition of the §6 item;
! NO hand-maintained reaction-type list except reading the code's own
! gtodn_jac_list for the value overlay):
!   bit 0 diagonal        row==col
!   bit 1 gas_chem        non-coag rxn, all-gas, col reactant & row compound
!   bit 2 surface_net     non-coag rxn touching a J or GRAIN species (folds
!                         grain-charging into surface -- design decision)
!   bit 3 coag_grain      COAGULATION_TYPE rxn, reactant-1 a GRAIN species
!   bit 4 ice_transport   COAGULATION_TYPE rxn, reactant-1 a J/ice species
!   bit 5 charge_electron cell involves e- of a coag rxn whose slot-6 is 'e-'
!   bit 6 grain_coupling  (row,col) in DECLARED_JAC_ROW/COL (live-divisor block)
!   bit 7 gtodn_recip     value overlay: r in gtodn_jac_list, bin k, compound
!                         row -> cell (row, INDGRAIN(k)/INDGRAIN_MINUS(k))
!
! Completeness gate (hard): every pattern cell must get mask/=0. A mask==0 cell
! is a classifier bug -- fail loudly with (row,col) and species names.
!
! Run from a prepared run directory (the fixture_jac_figure .in files). Exit 0 = PASS.
! ===========================================================================
program jacobian_provenance_dump
use iso_fortran_env, only: error_unit, output_unit
use global_variables
use ode_solver
use gasgrain
implicit none

integer :: N, p, col, row, r, s, d, k, i, a, bb, e_idx
integer :: nreac, ncomp
integer, dimension(MAX_COMPOUNDS) :: reac, comp
integer, allocatable :: mask(:)
logical :: is_coag, rxn_surface, has_eshed
integer :: bitcnt(0:7)
integer(8) :: totbits
integer :: mixed
integer :: pair(0:7,0:7)
integer :: rc_coag_grain, rc_ice_transport      ! reaction-level tallies
integer :: unclass_row, unclass_col, nzero
integer :: ua, ub
character(len=64) :: bitname(0:7)
character(len=8)  :: scls
character(len=2)  :: schg

bitname(0)='diagonal'      ; bitname(1)='gas_chem'
bitname(2)='surface_net'   ; bitname(3)='coag_grain'
bitname(4)='ice_transport' ; bitname(5)='charge_electron'
bitname(6)='grain_coupling'; bitname(7)='gtodn_recip'

call init_gasgrain()
N = nb_species

if (.not.allocated(JA_SYM) .or. NNZ_SYM <= 0) then
  write(error_unit,'(a)') 'FATAL: symbolic pattern not built (JA_SYM unallocated or NNZ_SYM<=0).'
  write(error_unit,'(a)') '  The fixture must set sparsity = symbolic and coagulation = 1.'
  call exit(2)
endif

e_idx = 0
do i=1,N
  if (trim(species_name(i)) == 'e-') e_idx = i
enddo

write(output_unit,'(a)') '=== jacobian_provenance_dump: init sizes ==='
write(output_unit,'(a,i0)')  '  nb_grains                 = ', nb_grains
write(output_unit,'(a,i0)')  '  nb_species (N)            = ', N
write(output_unit,'(a,i0)')  '  nb_reactions              = ', nb_reactions
write(output_unit,'(a,i0)')  '  nb_chemistry_reactions    = ', nb_chemistry_reactions
write(output_unit,'(a,i0)')  '  nb_coag_grain_reactions   = ', nb_coag_grain_reactions
write(output_unit,'(a,i0)')  '  nb_ice_transport_reactions= ', nb_ice_transport_reactions
write(output_unit,'(a,i0)')  '  nb_ice_bases              = ', nb_ice_bases
write(output_unit,'(a,i0)')  '  NB_DECLARED_JAC_DEPS      = ', NB_DECLARED_JAC_DEPS
write(output_unit,'(a,i0)')  '  gtodn_jac_n               = ', gtodn_jac_n
write(output_unit,'(a,i0)')  '  NNZ_SYM                   = ', NNZ_SYM
if (e_idx==0) write(output_unit,'(a)') '  (note: no e- species found; charge_electron bit will be empty)'

allocate(mask(NNZ_SYM)); mask = 0

! ---- bit 0: diagonal ----
do i=1,N
  p = locate(i,i)
  if (p>0) mask(p) = ibset(mask(p),0)
enddo

! ---- reactions: bits 1,2 (chemistry gas/surface), 3,4 (coag), 5 (charge e-) ----
rc_coag_grain = 0 ; rc_ice_transport = 0
do r=1,nb_reactions
  nreac=0 ; ncomp=0
  do s=1,MAX_REACTANTS
    if (valid(REACTION_COMPOUNDS_ID(s,r))) then
      nreac=nreac+1 ; reac(nreac)=REACTION_COMPOUNDS_ID(s,r)
    endif
  enddo
  do s=1,MAX_COMPOUNDS
    if (valid(REACTION_COMPOUNDS_ID(s,r))) then
      ncomp=ncomp+1 ; comp(ncomp)=REACTION_COMPOUNDS_ID(s,r)
    endif
  enddo
  if (nreac==0) cycle

  is_coag = (REACTION_TYPE(r)==COAGULATION_TYPE)

  if (is_coag) then
    ! reactant-1 decides grain vs ice (name-based; GRAIN* vs J*/K*)
    if (is_grain(reac(1))) then
      rc_coag_grain = rc_coag_grain + 1
      call set_bit_over_reaction(reac,nreac,comp,ncomp,3)
    else if (is_ice(reac(1))) then
      rc_ice_transport = rc_ice_transport + 1
      call set_bit_over_reaction(reac,nreac,comp,ncomp,4)
    else
      write(error_unit,'(a,i0,a,a)') 'FATAL: coag reaction ', r, &
        ' has reactant-1 neither GRAIN nor J/ice: ', trim(species_name(reac(1)))
      call exit(2)
    endif
    ! charge_electron overlay: slot-6 product == 'e-'
    has_eshed = .false.
    if (trim(REACTION_COMPOUNDS_NAMES(6,r)) == 'e-') has_eshed = .true.
    if (has_eshed .and. e_idx>0) then
      do a=1,nreac
        col = reac(a)
        do bb=1,ncomp
          row = comp(bb)
          if (row==e_idx .or. col==e_idx) then
            p = locate(row,col)
            if (p>0) mask(p) = ibset(mask(p),5)
          endif
        enddo
      enddo
    endif
  else
    ! chemistry: surface if the reaction touches any J or GRAIN species; else gas
    rxn_surface = .false.
    do s=1,ncomp
      if (is_ice(comp(s)) .or. is_grain(comp(s))) then
        rxn_surface = .true. ; exit
      endif
    enddo
    if (rxn_surface) then
      call set_bit_over_reaction(reac,nreac,comp,ncomp,2)
    else
      call set_bit_over_reaction(reac,nreac,comp,ncomp,1)
    endif
  endif
enddo

! ---- bit 6: declared non-reactant grain-coupling (live-divisor block) ----
do d=1,NB_DECLARED_JAC_DEPS
  p = locate(DECLARED_JAC_ROW(d), DECLARED_JAC_COL(d))
  if (p>0) mask(p) = ibset(mask(p),6)
enddo

! ---- bit 7: gtodn reciprocal value overlay (subset of bit 6) ----
do i=1,gtodn_jac_n
  r = gtodn_jac_list(i)
  k = GRAIN_RANK(r)
  if (k<1 .or. k>nb_grains) cycle
  do s=1,MAX_COMPOUNDS
    row = REACTION_COMPOUNDS_ID(s,r)
    if (.not.valid(row)) cycle
    if (INDGRAIN(k)>=1) then
      p = locate(row, INDGRAIN(k))
      if (p>0) mask(p) = ibset(mask(p),7)
    endif
    if (INDGRAIN_MINUS(k)>=1) then
      p = locate(row, INDGRAIN_MINUS(k))
      if (p>0) mask(p) = ibset(mask(p),7)
    endif
  enddo
enddo

! ====================== completeness gate (hard) ======================
nzero = 0 ; unclass_row=0 ; unclass_col=0
do col=1,N
  do p=IA_SYM(col), IA_SYM(col+1)-1
    if (mask(p)==0) then
      nzero = nzero + 1
      if (unclass_row==0) then
        unclass_row = JA_SYM(p) ; unclass_col = col
      endif
    endif
  enddo
enddo
if (nzero > 0) then
  write(error_unit,'(a,i0,a)') 'FATAL (completeness gate): ', nzero, ' pattern cell(s) with mask==0.'
  write(error_unit,'(a,i0,a,i0,a)') '  first unclassifiable cell: (row=', unclass_row, &
    ', col=', unclass_col, ')'
  write(error_unit,'(4a)') '    row species = ', trim(species_name(unclass_row)), &
    '   col species = ', trim(species_name(unclass_col))
  call exit(3)
endif

! consistency: every bit-7 cell must also be bit-6
do p=1,NNZ_SYM
  if (btest(mask(p),7) .and. .not.btest(mask(p),6)) then
    write(error_unit,'(a)') 'FATAL: a gtodn_recip (bit7) cell is not in the declared block (bit6).'
    call exit(3)
  endif
enddo

! ====================== write TSVs ======================
open(unit=71, file='jac_pattern.tsv', status='replace', action='write')
write(71,'(a)') 'row'//char(9)//'col'//char(9)//'mech_mask'
do col=1,N
  do p=IA_SYM(col), IA_SYM(col+1)-1
    row = JA_SYM(p)
    write(71,'(i0,a,i0,a,i0)') row, char(9), col, char(9), mask(p)
  enddo
enddo
close(71)

open(unit=72, file='jac_species.tsv', status='replace', action='write')
write(72,'(a)') 'index'//char(9)//'name'//char(9)//'class'//char(9)//'bin'//char(9)//'charge'
do i=1,N
  ! class
  if (is_grain(i)) then
    scls = 'grain'
  else if (is_ice(i)) then
    scls = 'ice'
  else
    scls = 'gas'
  endif
  ! charge (grains only: '-' if name carries a trailing minus, else '0')
  schg = ''
  if (is_grain(i)) then
    if (index(trim(species_name(i)),'-') > 0) then
      schg = '-'
    else
      schg = '0'
    endif
  endif
  write(72,'(i0,a,a,a,a,a,i0,a,a)') i, char(9), trim(species_name(i)), char(9), &
    trim(scls), char(9), species_bin(i), char(9), trim(schg)
enddo
close(72)

! ====================== gate report ======================
bitcnt = 0 ; totbits = 0_8
do p=1,NNZ_SYM
  do i=0,7
    if (btest(mask(p),i)) bitcnt(i) = bitcnt(i) + 1
  enddo
  totbits = totbits + popcnt(mask(p))
enddo

mixed = 0
pair = 0
do p=1,NNZ_SYM
  if (popcnt(mask(p)) >= 2) mixed = mixed + 1
  do a=0,7
    if (.not.btest(mask(p),a)) cycle
    do bb=a+1,7
      if (btest(mask(p),bb)) pair(a,bb) = pair(a,bb) + 1
    enddo
  enddo
enddo

write(output_unit,'(a)') ''
write(output_unit,'(a)') '=== COMPLETENESS GATE: PASS (no mask==0 cell) ==='
write(output_unit,'(a,i0)') 'NNZ_SYM                    = ', NNZ_SYM
write(output_unit,'(a,i0)') 'total set bits (sum popcnt) = ', totbits
write(output_unit,'(a)') ''
write(output_unit,'(a)') 'per-mechanism CELL counts (bit : name : cells):'
do i=0,7
  write(output_unit,'(a,i1,a,a16,a,i0)') '  bit ', i, '  ', bitname(i), ' : ', bitcnt(i)
enddo
write(output_unit,'(a)') ''
write(output_unit,'(a,i0,a,f5.1,a)') 'mixed cells (>=2 bits)     = ', mixed, &
  '   (', 100.0*real(mixed)/real(NNZ_SYM), ' % of NNZ)'
write(output_unit,'(a)') 'overlapping bit-pairs (nonzero, cells):'
do a=0,7
  do bb=a+1,7
    if (pair(a,bb) > 0) then
      write(output_unit,'(a,a16,a,a16,a,i0)') '  ', bitname(a), ' & ', bitname(bb), ' : ', pair(a,bb)
    endif
  enddo
enddo
write(output_unit,'(a)') ''
write(output_unit,'(a)') 'reaction-level cross-check (tags vs code counters):'
write(output_unit,'(a,i0,a,i0)') '  coag_grain reactions    tagged=', rc_coag_grain, &
  '  code nb_coag_grain_reactions=', nb_coag_grain_reactions
write(output_unit,'(a,i0,a,i0)') '  ice_transport reactions tagged=', rc_ice_transport, &
  '  code nb_ice_transport_reactions=', nb_ice_transport_reactions
if (rc_coag_grain /= nb_coag_grain_reactions .or. rc_ice_transport /= nb_ice_transport_reactions) then
  write(error_unit,'(a)') 'FATAL: reaction-level tag counts disagree with the code counters.'
  call exit(3)
endif
write(output_unit,'(a)') ''
write(output_unit,'(a,i0)') 'cross-check band: coag_grain cells    (expect low hundreds)      = ', bitcnt(3)
write(output_unit,'(a,i0)') 'cross-check band: ice_transport cells (expect high 100s-low 1000s)= ', bitcnt(4)
write(output_unit,'(a,i0)') 'cross-check band: grain_coupling cells(expect 2544)              = ', bitcnt(6)

write(output_unit,'(a)') ''
write(output_unit,'(a)') 'wrote jac_pattern.tsv and jac_species.tsv'
write(output_unit,'(a)') 'DONE (gate PASS).'

contains

  logical function valid(idx)
    integer, intent(in) :: idx
    valid = (idx>=1 .and. idx<=nb_species)
  end function valid

  ! binary search for row in column col of CSC pattern; 0 if absent
  integer function locate(irow, icol) result(pos)
    integer, intent(in) :: irow, icol
    integer :: lo, hi, mid
    pos = 0
    if (icol<1 .or. icol>nb_species) return
    lo = IA_SYM(icol) ; hi = IA_SYM(icol+1)-1
    do while (lo<=hi)
      mid = (lo+hi)/2
      if (JA_SYM(mid)==irow) then
        pos = mid ; return
      else if (JA_SYM(mid)<irow) then
        lo = mid+1
      else
        hi = mid-1
      endif
    enddo
  end function locate

  subroutine set_bit_over_reaction(rlist, nr, clist, nc, ibit)
    integer, intent(in) :: rlist(:), clist(:), nr, nc, ibit
    integer :: aa, bbi, pp
    do aa=1,nr
      do bbi=1,nc
        pp = locate(clist(bbi), rlist(aa))
        if (pp>0) mask(pp) = ibset(mask(pp), ibit)
      enddo
    enddo
  end subroutine set_bit_over_reaction

  logical function is_grain(idx)
    integer, intent(in) :: idx
    character(len=len(species_name(1))) :: nm
    is_grain = .false.
    if (idx<1 .or. idx>nb_species) return
    nm = species_name(idx)
    if (len_trim(nm)>=5) then
      if (nm(1:5)=='GRAIN') is_grain = .true.
    endif
  end function is_grain

  logical function is_ice(idx)
    integer, intent(in) :: idx
    character(len=len(species_name(1))) :: nm
    is_ice = .false.
    if (idx<1 .or. idx>nb_species) return
    nm = species_name(idx)
    if (nm(1:1)=='J' .or. nm(1:1)=='K') is_ice = .true.
  end function is_ice


  ! grain bin from grain_col_bin (grains) or parsed 2-digit index (ice Jkk..), else 0
  integer function species_bin(idx)
    integer, intent(in) :: idx
    character(len=len(species_name(1))) :: nm
    integer :: ios, kk
    species_bin = 0
    if (idx<1 .or. idx>nb_species) return
    if (allocated(grain_col_bin)) then
      if (grain_col_bin(idx) > 0) then
        species_bin = grain_col_bin(idx) ; return
      endif
    endif
    nm = species_name(idx)
    if (is_ice(idx) .and. len_trim(nm)>=3) then
      read(nm(2:3),'(i2)',iostat=ios) kk
      if (ios==0) species_bin = kk
    endif
  end function species_bin


end program jacobian_provenance_dump
