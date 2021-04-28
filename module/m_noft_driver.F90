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
 use m_optocc

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

subroutine run_noft(INOF_in,Ista_in,NBF_occ_in,Nfrozen_in,Npairs_in,&
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in,imethocc,Vnn,hCORE,ERI_J,ERI_K)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::INOF_in,Ista_in,imethocc
 integer,intent(in)::NBF_occ_in,Nfrozen_in,Npairs_in,Ncoupled_in
 integer,intent(in)::Nbeta_elect_in,Nalpha_elect_in
 double precision,intent(in)::Vnn
!arrays
 double precision,dimension(:),intent(inout)::hCORE,ERI_J,ERI_K
!Local variables ------------------------------
!scalars
 type(rdm_t),target::RDMd
!arrays

 call rdm_init(RDMd,INOF_in,Ista_in,NBF_occ_in,Nfrozen_in,Npairs_in,&
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in)

 !TODO
 call opt_occ(imethocc,RDMd,Vnn,hCORE,ERI_J,ERI_K)


 call RDMd%free() 

end subroutine run_noft
!!***

end module m_noft_driver
!!***
