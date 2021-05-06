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
!! GAMMAs=Independent variables used in the unconstrained optimization as cos^2 (gamma) + sin^2 (gamma) = 1
!!
!! PARENTS
!!
!! CHILDREN
!!   m_optocc
!!
!! SOURCE
module m_noft_driver

 use m_rdmd
 use m_integd
 use m_elag
 use m_optocc
 use m_optorb

 implicit none

!! private :: 
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
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in,imethocc,itermax,iprintdmn,&
&  tolE,Vnn,NO_COEF,mo_ints)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::INOF_in,Ista_in,imethocc,itermax,iprintdmn
 integer,intent(in)::NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,Ncoupled_in
 integer,intent(in)::Nbeta_elect_in,Nalpha_elect_in
 double precision,intent(in)::Vnn,tolE
 external::mo_ints
!arrays
 double precision,dimension(NBF_tot_in,NBF_tot_in),intent(inout)::NO_COEF
!Local variables ------------------------------
!scalars
 logical::ekt
 integer::iorb,iorb1,iter
 double precision::Energy,Energy_old
 type(rdm_t),target::RDMd
 type(integ_t),target::INTEGd
!arrays
 character(len=10)::coef_file
!************************************************************************

 call rdm_init(RDMd,INOF_in,Ista_in,NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in)
 call integ_init(INTEGd,RDMd%NBF_tot,RDMd%NBF_occ)

 ! Occ optimization using guess orbs. (HF, CORE, etc).
 write(*,'(a)') ' '
 iter=0
 call mo_ints(NO_COEF,INTEGd%hCORE,INTEGd%ERImol)
 call INTEGd%htohpp(RDMd%NBF_occ)
 call INTEGd%eritoeriJK(RDMd%NBF_occ)
 call opt_occ(iter,imethocc,RDMd,Vnn,Energy,INTEGd%hCOREpp,INTEGd%ERI_J,INTEGd%ERI_K)
 Energy_old=Energy

 ! Orb. and occ. optimization
 do
  ! Orb. optimization
  call opt_orb(RDMd,INTEGd,Vnn,Energy,NO_COEF,mo_ints)

  ! Occ. optimization
  call opt_occ(iter,imethocc,RDMd,Vnn,Energy,INTEGd%hCOREpp,INTEGd%ERI_J,INTEGd%ERI_K)

  ! Check convergence
  if(abs(Energy-Energy_old)<tolE) then
   Energy_old=Energy
   exit
  endif
  Energy_old=Energy
  
  ! Check maximum number of iterations
  if(iter>itermax) exit

 enddo

 if(iprintdmn==1) call RDMd%print_dmn(RDMd%DM2_J,RDMd%DM2_K) 

 ! Free all allocated INTEGd arrays
 call INTEGd%free()

 ! Print final diagonalized INTEGd%Lambdas values
 call diag_ekt(RDMd,NO_COEF)

 ! Print final EKT values
 call diag_ekt(RDMd,NO_COEF,ekt)

 ! Print final occ. numbers
 write(*,'(a)') ' '
 RDMd%OCC(:)=2.0d0*RDMd%OCC(:)
 write(*,'(a,f10.5,a)') 'Total occ ',sum(RDMd%OCC(:)),' final occ. numbers '
 iorb1=RDMd%NBF_occ-(RDMd%NBF_occ/10)*10
 do iorb=1,(RDMd%NBF_occ/10)*10,10
  write(*,'(f12.6,9f11.6)') RDMd%OCC(iorb:iorb+9)
 enddo
 iorb1=(RDMd%NBF_occ/10)*10+1 
 write(*,'(f12.6,*(f11.6))') RDMd%OCC(iorb1:) 
 write(*,'(a)') ' '

 ! Print final nat. orb. coef.
 coef_file='NO_COEF'
 call RDMd%print_orbs(NO_COEF,coef_file)

 ! Free all allocated RDMd arrays
 call RDMd%free() 

end subroutine run_noft
!!***

end module m_noft_driver
!!***
