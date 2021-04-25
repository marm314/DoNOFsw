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
 use m_E_grad_occ

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
!!  Run optimization w.r.t occ numbers 
!!
!! INPUTS
!! HighSpin_in=Logical variable to decide what spin-uncompensated version to use (default=False i.e. use mixture of states)
!! MSpin_in=Integer variable to define the MS (default=0 i.e. spin-compensated)
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

subroutine run_noft(HighSpin_in,MSpin_in,INOF_in,Ista_in,NBF_occ_in,Nfrozen_in,Npairs_in,&
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in,Nvars,hCORE,ERI_J,ERI_K,GAMMAs_in,RO,CJ12,CK12,DR,DCJ12r,DCK12r)
!Arguments ------------------------------------
!scalars
 logical,intent(in)::HighSpin_in
 integer,intent(in)::INOF_in,MSpin_in,Ista_in
 integer,intent(in)::NBF_occ_in,Nfrozen_in,Npairs_in,Ncoupled_in
 integer,intent(in)::Nbeta_elect_in,Nalpha_elect_in
!arrays
 double precision,dimension(:),intent(inout)::hCORE,ERI_J,ERI_K
 
 !!!!!!!!!!!!!!!
 !! TO REMOVE !!
 !!!!!!!!!!!!!!!
 integer,intent(in)::Nvars
 double precision,dimension(Nvars),intent(inout)::GAMMAs_in
 double precision,dimension(NBF_occ_in),intent(inout)::RO
 double precision,dimension(NBF_occ_in*NBF_occ_in),intent(inout)::CJ12,CK12
 double precision,dimension(NBF_occ_in*Nvars),intent(inout)::DR
 double precision,dimension(NBF_occ_in*NBF_occ_in*Nvars),intent(inout)::DCJ12r,DCK12r
 integer::ind
 !!!!!!!!!!!!!!!

!Local variables ------------------------------
!scalars
 double precision::Energy
 double precision,dimension(:),allocatable::GAMMAs
 type(rdm_t),target::RDMd
!arrays

 call rdm_init(RDMd,HighSpin_in,MSpin_in,INOF_in,Ista_in,NBF_occ_in,Nfrozen_in,Npairs_in,&
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in)
 allocate(GAMMAs(RDMd%Ngammas))

 do ind=1,NBF_occ_in
  RDMd%occ(ind)=RO(ind)
 enddo
 do ind=1,NBF_occ_in*NBF_occ_in
  RDMd%DM2_J(ind)=CJ12(ind)
  RDMd%DM2_K(ind)=CK12(ind)
 enddo
 do ind=1,NBF_occ_in*Nvars
  RDMd%Docc_gamma(ind)=DR(ind)
 enddo
 do ind=1,NBF_occ_in*NBF_occ_in*Nvars
  RDMd%DDM2_gamma_J(ind)=DCJ12r(ind)
  RDMd%DDM2_gamma_K(ind)=DCK12r(ind)
 enddo

 GAMMAs=GAMMAs_in
 call calc_E_occ(RDMd,GAMMAs,Energy,hCORE,ERI_J,ERI_K)
 write(*,*) Energy

 do ind=1,NBF_occ_in
  RO(ind)=RDMd%occ(ind)
 enddo
 do ind=1,NBF_occ_in*NBF_occ_in
  CJ12(ind)=RDMd%DM2_J(ind)     
  CK12(ind)=RDMd%DM2_K(ind)
 enddo 
 do ind=1,NBF_occ_in*Nvars
  DR(ind)=RDMd%Docc_gamma(ind)
 enddo
 do ind=1,NBF_occ_in*NBF_occ_in*Nvars
  DCJ12r(ind)=RDMd%DDM2_gamma_J(ind)
  DCK12r(ind)=RDMd%DDM2_gamma_K(ind)
 enddo

 deallocate(GAMMAs)
 call RDMd%free() 

end subroutine run_noft
!!***

end module m_noft_driver
!!***
