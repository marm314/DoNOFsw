!!****m* DoNOF/m_noft_driver
!! NAME
!!  m_noft_driver
!!
!! FUNCTION
!!  Module prepared to perform all procedures required for occ. and orbital optmization  
!!
!!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  Nfrozen |            Npairs_p_sing     |              Nvirt                               = NBF
!!  Nfrozen |         Npairs + Nsingleocc  |     Ncoupled*Npairs                   + Nempty   = NBF
!!                           | Nsingleocc  |   NBF_occ - Npairs_p_sing - Nfrozen   | Nempty   = NBF
!!                           Nbeta         Nalpha                                  NBF_occ
!!- - - - - - - - - - - - - - -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!!
!! PARENTS
!!
!! CHILDREN
!!   m_optocc
!!   m_optorb
!!
!! SOURCE
module m_noft_driver

 use m_rdmd
 use m_integd
 use m_elag
 use m_optocc
 use m_optorb

 implicit none

 private :: read_restart,echo_input
!!***

 public :: run_noft
!!***

contains

!!***
!!****f* DoNOF/run_noft
!! NAME
!! run_noft
!!
!! FUNCTION
!!  Run NOFT procedures 
!!
!! INPUTS
!! INOF_in=PNOFi functional to use
!! Ista_in=Use PNOF7 (Ista_in=0) or PNOF7s (Ista_in=1)
!! NBF_occ_in=Number of orbitals that are occupied
!! Nfrozen_in=Number of frozen orbitals that remain with occ=2.0 
!! Npairs_in=Number of electron pairs
!! Ncoupled_in=Number of coupled orbitals per electron pair MINUS ONE
!! Nbeta_elect_in=Number of beta electrons (N/2 for spin compensated systems)
!! Nalpha_elect_in=Number of beta electrons (N/2 for spin compensated systems)
!!
!! OUTPUT
!! occ=Occupancies of the frozen + active orbitals
!! DM2_J=DM2 elements that use J integrals 
!! DM2_K=DM2 elements that use K integrals 
!! DDM2_gamma_J=Derivative of the DM2 elements w.r.t. gamma that use J integrals 
!! DDM2_gamma_K=Derivative of the DM2 elements w.r.t. gamma that use K integrals
!! Docc_gamma=Derivative of the occupancies w.r.t. gamma
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine run_noft(INOF_in,Ista_in,NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,&
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in,imethocc,imethorb,itermax,iprintdmn,iprintints,&
&  itolLambda,ndiis,tolE_in,Vnn,NO_COEF,Overlap_in,mo_ints,restart,ireadGAMMAS,ireadCOEF,&
&  ireadFdiag)
!Arguments ------------------------------------
!scalars
 integer,optional,intent(in)::ireadGAMMAS,ireadCOEF,ireadFdiag
 logical,optional,intent(in)::restart
 integer,intent(in)::INOF_in,Ista_in,imethocc,imethorb,itermax,iprintdmn,iprintints,itolLambda,ndiis
 integer,intent(in)::NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,Ncoupled_in
 integer,intent(in)::Nbeta_elect_in,Nalpha_elect_in
 double precision,intent(in)::Vnn,tolE_in
 external::mo_ints
!arrays
 double precision,dimension(NBF_tot_in,NBF_tot_in),intent(in)::Overlap_in
 double precision,dimension(NBF_tot_in,NBF_tot_in),intent(inout)::NO_COEF
!Local variables ------------------------------
!scalars
 logical::ekt,diagLpL,restart_param
 integer::iorb,iter
 double precision::Energy,Energy_old,Vee,hONEbody
 type(rdm_t),target::RDMd
 type(integ_t),target::INTEGd
 type(elag_t),target::ELAGd
!arrays
 character(len=10)::coef_file
!************************************************************************

 diagLpL=.true.; restart_param=.false.;

 ! Write Header
 write(*,'(a)') ' '
 write(*,'(a)') ' -------------------------------------------'
 write(*,'(a)') ' Entering RUN-NOF module for NOFT calcs.'
 write(*,'(a)') ' '
 write(*,'(a)') ' Developed by: Dr. M. Rodriguez-Mayorga '
 write(*,'(a)') ' '
 write(*,'(a)') '  First version: VU Amsterdam 2021 '
 write(*,'(a)') ' '
 write(*,'(a)') ' -------------------------------------------'
 write(*,'(a)') ' '

 ! Print user defined parameters used in this run
 if(present(restart)) then
  if(present(ireadGAMMAS).and.present(ireadCOEF).and.present(ireadFdiag)) then
   restart_param=.true.
   call echo_input(INOF_in,Ista_in,NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,&
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in,imethocc,imethorb,itermax,iprintdmn,&
&  iprintints,itolLambda,ndiis,tolE_in,restart=restart,ireadGAMMAS=ireadGAMMAS,&
&  ireadCOEF=ireadCOEF,ireadFdiag=ireadFdiag)
  else
   write(*,'(a)') 'Warning! Asking for restart but the restart parameters are unspecified (not restarting).' 
   restart_param=.false.
   call echo_input(INOF_in,Ista_in,NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,&
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in,imethocc,imethorb,itermax,iprintdmn,&
&  iprintints,itolLambda,ndiis,tolE_in)
  endif
 else 
  restart_param=.false.
  call echo_input(INOF_in,Ista_in,NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,&
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in,imethocc,imethorb,itermax,iprintdmn,&
&  iprintints,itolLambda,ndiis,tolE_in)
 endif

 ! Initialize RDMd, INTEGd, and ELAGd objects.
 call rdm_init(RDMd,INOF_in,Ista_in,NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in)
 call integ_init(INTEGd,RDMd%NBF_tot,RDMd%NBF_occ,Overlap_in)
 call elag_init(ELAGd,RDMd%NBF_tot,diagLpL,itolLambda,ndiis,imethorb,tolE_in)

 ! Check for the presence of restart files. If they are available read them if demanded 
 if(restart_param) then
  write(*,'(a)') ' '
  call read_restart(RDMd,ELAGd,NO_COEF,ireadGAMMAS,ireadCOEF,ireadFdiag)
  write(*,'(a)') ' '
 endif

 ! Occ optimization using guess orbs. (HF, CORE, etc).
 write(*,'(a)') ' '
 iter=-1;
 call mo_ints(NO_COEF,INTEGd%hCORE,INTEGd%ERImol)
 call INTEGd%eritoeriJK(RDMd%NBF_occ)
 call opt_occ(iter,imethocc,RDMd,Vnn,Energy,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K) ! Also iter=iter+1
 Energy_old=Energy

 ! Orb. and occ. optimization
 do
  ! Orb. optimization
  call ELAGd%clean_diis()
  call opt_orb(iter,imethorb,ELAGd,RDMd,INTEGd,Vnn,Energy,NO_COEF,mo_ints)
  if(imethorb==1) then ! For F diag method, print F_pp elements after each global iteration
   call ELAGd%print_Fdiag(RDMd%NBF_tot)
  endif
  call RDMd%print_orbs_bin(NO_COEF)

  ! Occ. optimization
  call opt_occ(iter,imethocc,RDMd,Vnn,Energy,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K) ! Also iter=iter+1
  call RDMd%print_gammas()

  ! Check convergence
  if(abs(Energy-Energy_old)<ELAGd%tolE) then
   Energy_old=Energy
   exit
  endif
  Energy_old=Energy
  
  ! Check maximum number of iterations
  if(iter>itermax) exit

 enddo

 if(iprintdmn==1) call RDMd%print_dmn(RDMd%DM2_J,RDMd%DM2_K) 

 ! Print final diagonalized INTEGd%Lambdas values
 call ELAGd%diag_lag(RDMd,INTEGd,NO_COEF)

 ! Print final Extended Koopmans' Theorem (EKT) values
 if(RDMd%Nsingleocc==0) call ELAGd%diag_lag(RDMd,INTEGd,NO_COEF,ekt=ekt)

 ! Print final occ. numbers
 write(*,'(a)') ' '
 RDMd%occ(:)=2.0d0*RDMd%occ(:)
 write(*,'(a,f10.5,a)') 'Total occ ',sum(RDMd%occ(:)),' optimized occ. numbers '
 do iorb=1,(RDMd%NBF_occ/10)*10,10
  write(*,'(f12.6,9f11.6)') RDMd%occ(iorb:iorb+9)
 enddo
 iorb=(RDMd%NBF_occ/10)*10+1 
 write(*,'(f12.6,*(f11.6))') RDMd%occ(iorb:) 
 write(*,'(a)') ' '

 ! Print final nat. orb. coef.
 coef_file='NO_COEF'
 call RDMd%print_orbs(NO_COEF,coef_file)
 call RDMd%print_orbs_bin(NO_COEF)
 
 ! Print final Energy and its components
 hONEbody=0.0d0
 do iorb=1,RDMd%NBF_occ
  hONEbody=hONEbody+RDMd%occ(iorb)*INTEGd%hCORE(iorb,iorb)
 enddo
 Vee=Energy-hONEbody
 write(*,'(a)') ' '
 write(*,'(a,f15.6,a,i6,a)') 'Final NOF energy= ',Energy+Vnn,' a.u. after ',iter,' global iter.'
 write(*,'(a,f15.6,a)') 'hCORE           = ',hONEbody,' a.u.'
 write(*,'(a,f15.6,a)') 'Vee             = ',Vee,' a.u.'
 write(*,'(a,f15.6,a)') 'Vnn             = ',Vnn,' a.u.'
 write(*,'(a)') ' '

 ! Print the hCORE and ERImol ints in the latest NO_COEF basis
 if(iprintints==1) then
  call INTEGd%print_int(RDMd%NBF_tot)
 endif

 ! Free all allocated RDMd, INTEGd, and ELAGd arrays
 call ELAGd%free() 
 call INTEGd%free()
 call RDMd%free() 

 ! Write Footer
 write(*,'(a)') ' '
 write(*,'(a)') ' -------------------------------------------'
 write(*,'(a)') ' '
 write(*,'(a)') ' Normal termination of RUN-NOF module.'
 write(*,'(a)') ' '
 write(*,'(a)') '   |                   | '
 write(*,'(a)') '   |                   | '
 write(*,'(a)') '   |                   | '
 write(*,'(a)') '   |        <^>        | '
 write(*,'(a)') '   ||===I||(-@-)||I===|| '
 write(*,'(a)') '   |        \_/        | '
 write(*,'(a)') '   |                   | '
 write(*,'(a)') '   |                   | '
 write(*,'(a)') '   |                   | '
 write(*,'(a)') '   |                   | '
 write(*,'(a)') '   |                   | '
 write(*,'(a)') ' '
 write(*,'(a)') ' "Your feeble skills are no match for the '
 write(*,'(a)') ' power of the dark side." Emperor Palpatine '
 write(*,'(a)') ' '
 write(*,'(a)') ' -------------------------------------------'
 write(*,'(a)') ' '

end subroutine run_noft
!!***

!!***
!!****f* DoNOF/echo_input
!! NAME
!! echo_input
!!
!! FUNCTION
!!  Echo all parameters employed in this run of this module
!!
!! INPUTS
!! INOF_in=PNOFi functional to use
!! Ista_in=Use PNOF7 (Ista_in=0) or PNOF7s (Ista_in=1)
!! NBF_occ_in=Number of orbitals that are occupied
!! Nfrozen_in=Number of frozen orbitals that remain with occ=2.0 
!! Npairs_in=Number of electron pairs
!! Ncoupled_in=Number of coupled orbitals per electron pair MINUS ONE
!! Nbeta_elect_in=Number of beta electrons (N/2 for spin compensated systems)
!! Nalpha_elect_in=Number of beta electrons (N/2 for spin compensated systems)
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine echo_input(INOF_in,Ista_in,NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,&
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in,imethocc,imethorb,itermax,iprintdmn,&
&  iprintints,itolLambda,ndiis,tolE_in,restart,ireadGAMMAS,ireadCOEF,ireadFdiag)
!Arguments ------------------------------------
!scalars
 logical,optional,intent(in)::restart
 integer,optional,intent(in)::ireadGAMMAS,ireadCOEF,ireadFdiag
 integer,intent(in)::INOF_in,Ista_in,imethocc,imethorb,itermax,iprintdmn,iprintints,itolLambda,ndiis
 integer,intent(in)::NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,Ncoupled_in
 integer,intent(in)::Nbeta_elect_in,Nalpha_elect_in
 double precision,intent(in)::tolE_in
!arrays
!Local variables ------------------------------
!scalars
!arrays
!************************************************************************
 write(*,'(a)') ' '
 write(*,'(a,i12)') ' NOF approximation in use          ',INOF_in
 if(INOF_in==7) then
  if(Ista_in==1) then
   write(*,'(a,i12)') ' PNOF7s version selected Istat     ',Ista_in
  else
   write(*,'(a,i12)') ' PNOF7  version selected Istat     ',Ista_in
  endif
 endif
 write(*,'(a,i12)') ' Numb. of basis functions          ',NBF_tot_in
 write(*,'(a,i12)') ' Numb. of occ orbitals             ',NBF_occ_in
 write(*,'(a,i12)') ' Numb. of frozen orbs (occ=2)      ',NBF_occ_in
 write(*,'(a,i12)') ' Numb. of active e- pairs          ',Npairs_in
 write(*,'(a,i12)') ' Numb. of "virtual" coupled orbs   ',Ncoupled_in
 write(*,'(a,i12)') ' Numb. of singly occupied orbs     ',Nalpha_elect_in-Nbeta_elect_in
 if(imethocc==1) then
  write(*,'(a,i12)') ' CG method used in occ opt.        ',imethocc
 else
  write(*,'(a,i12)') ' L-BFGS method used in occ opt.    ',imethocc
 endif
 if(imethorb==1) then
  write(*,'(a,i12)') ' F_diag method used in orb opt.    ',imethorb
  write(*,'(a,e10.3)') ' Tolerance Lambda convergence        ',1.0d1**(-itolLambda)
  write(*,'(a,i12)') ' Numb. of iter used in DIIS        ',ndiis
 else
  write(*,'(a,i12)') ' Newton method used in orb opt.    ',imethorb
 endif
 write(*,'(a,i11)') ' Max. number of global iterations   ',itermax
 write(*,'(a,e10.3)') ' Tolerance Energy convergence        ',tolE_in
 write(*,'(a,i12)') ' Print optimal 1,2-RDMs (true=1)   ',iprintdmn
 write(*,'(a,i12)') ' Print last hCORE and ERImol ints  ',iprintints
 ! Check for the presence of restart files. If they are available, read them if required (default=not to read)
 if(present(restart)) then
  write(*,'(a,i12)') ' Restart reading GAMMAs (true=1)   ',ireadGAMMAS
  write(*,'(a,i12)') ' Restart reading COEFs  (true=1)   ',ireadCOEF
  if(imethorb==1) then
   write(*,'(a,i12)') ' Restart reading F_pp   (true=1)   ',ireadFdiag
  endif
 endif
 write(*,'(a)') ' '
 
end subroutine echo_input
!!***

!!***
!!****f* DoNOF/read_restart
!! NAME
!!  read_restart
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine read_restart(RDMd,ELAGd,NO_COEF,ireadGAMMAS,ireadCOEF,ireadFdiag)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::ireadGAMMAS,ireadCOEF,ireadFdiag
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
!arrays
 double precision,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::NO_COEF
!Local variables ------------------------------
!scalars
 integer::iunit,istat,intvar,intvar1,icount
 double precision::doubvar
 double precision,allocatable,dimension(:)::GAMMAS_in
 double precision,allocatable,dimension(:,:)::NO_COEF_in
!arrays
!************************************************************************

 ! Read NO_COEF for guess
 allocate(NO_COEF_in(RDMd%NBF_tot,RDMd%NBF_tot))
 open(unit=iunit,form='unformatted',file='NO_COEF_BIN',iostat=istat,status='old')
 icount=0
 if(istat==0.and.ireadCOEF==1) then
  do
   read(iunit,iostat=istat) intvar,intvar1,doubvar
   if(istat/=0) then
    exit
   endif
   if(((intvar/=0).and.(intvar1/=0)).and.intvar*intvar1<=RDMd%NBF_tot*RDMd%NBF_tot) then
    NO_COEF_in(intvar,intvar1)=doubvar
    icount=icount+1
   else
    exit
   endif
  enddo
 endif
 if(icount==RDMd%NBF_tot*RDMd%NBF_tot) then
  NO_COEF(:,:)=NO_COEF_in(:,:)
  write(*,'(a)') 'NO coefs. read from checkpoint file'
 endif
 close(iunit)
 deallocate(NO_COEF_in)

 ! Read GAMMAs indep. parameters used to compute occs.
 allocate(GAMMAS_in(RDMd%Ngammas))
 open(unit=iunit,form='unformatted',file='GAMMAS',iostat=istat,status='old')
 icount=0
 if(istat==0.and.ireadGAMMAS==1) then
  do 
   read(iunit,iostat=istat) intvar,doubvar
   if(istat/=0) then
    exit
   endif
   if((intvar/=0).and.intvar<=RDMd%Ngammas) then
    GAMMAs_in(intvar)=doubvar
    icount=icount+1
   else
    exit
   endif
  enddo
 endif
 if(icount==RDMd%Ngammas) then
  RDMd%GAMMAs_old(:)=GAMMAS_in(:)
  RDMd%GAMMAs_nread=.false.
  write(*,'(a)') 'GAMMAs (indep. variables) read from checkpoint file'
 endif
 close(iunit)
 deallocate(GAMMAS_in)

 ! Read diag. part of F matrix
 if(ELAGd%imethod==1) then
  open(unit=iunit,form='unformatted',file='F_DIAG',iostat=istat,status='old')
  icount=0
  if(istat==0.and.ireadFdiag==1) then
   do 
    read(iunit,iostat=istat) intvar,doubvar
    if(istat/=0) then
     exit
    endif
    if((intvar/=0).and.intvar<=RDMd%NBF_tot) then
     ELAGd%F_diag(intvar)=doubvar
     icount=icount+1
    else
     exit
    endif
   enddo
  endif
  if(icount==RDMd%NBF_tot) then
   ELAGd%diagLpL=.false.
   write(*,'(a)') 'F_pp elements read from checkpoint file'
  endif
  close(iunit)
 endif

end subroutine read_restart
!!***

end module m_noft_driver
!!***
