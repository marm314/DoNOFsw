!!****m* DoNOF/m_rdmd
!! NAME
!! basic NOFT variables
!!
!! FUNCTION
!! This module contains definitions of common variables used by NOFT module.
!!
!! COPYRIGHT
!! This file is distributed under the terms of the
!! GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! SOURCE

module m_rdmd 

 implicit none

!!****t* m_rdmd/rdm_t
!! NAME
!! rdm_t
!!
!! FUNCTION
!! Datatype storing noft quantities and arrays needed
!!
!! SOURCE

 type,public :: rdm_t

  logical::HighSpin=.false.      ! Decide if it is a high-spin calc. or mixture of states
  integer::INOF=7                ! Functional to use (5-> PNOF5, 7-> PNOF7)
  integer::MSpin=0               ! MS spin value
  integer::Ista=1                ! Use PNOF7s version
  integer::Nfrozen               ! Number of frozen orbitals in the NOFT calc.
  integer::Nbeta_elect           ! Number of orbitals containing beta electrons
  integer::Nalpha_elect          ! Number of orbitals containing alpha electrons
  integer::Nsingleocc=0          ! Number of singly occ orbitals
  integer::NBF_occ               ! Number of frozen plus active orbitals
  integer::NBF_ldiag             ! Size of the arrays that contain J and K integrals
  integer::Ncoupled              ! Number of 'virtual' coupled orbital per 'occupied' orbital  
  integer::Npairs                ! Number of electron pairs
  integer::Npairs_p_sing         ! Number of electron pairs plus number of singly occ orbitals
  integer::Ngammas               ! Number of gammas (independet variables used in occ optimization procedure)
  double precision::Sums
! arrays 
  double precision,allocatable,dimension(:)::GAMMAs
  double precision,allocatable,dimension(:)::OCC
  double precision,allocatable,dimension(:)::DM2_J,DM2_K
  double precision,allocatable,dimension(:)::Docc_gamma
  double precision,allocatable,dimension(:)::DDM2_gamma_J,DDM2_gamma_K

 contains 
   procedure :: free => rdm_free
   ! Destructor.

 end type rdm_t

 public :: rdm_init     ! Main creation method.
!!***

CONTAINS  !==============================================================================

!!***
!!****f* DoNOF/rdm_init
!! NAME
!! rdm_init
!!
!! FUNCTION
!!  Initialize the data type rdm_t 
!!
!! INPUTS
!! HighSpin=Logical variable to decide what spin-uncompensated version to use (default=False i.e. use mixture of states)
!! MSpin=Integer variable to define the MS (default=0 i.e. spin-compensated)
!! INOF=PNOFi functional to use
!! Ista=Use PNOF7 (Ista=0) or PNOF7s (Ista=1)
!! NBF_occ=Number of orbitals that are occupied
!! Nfrozen=Number of frozen orbitals that remain with occ=2.0 
!! Npairs=Number of electron pairs
!! Ncoupled=Number of coupled orbitals per electron pair (it is then used as Ncoupled-1 inside this module, as the number
!!             of coupled 'virtual' orbitals to a 'initially occupied' (HF) orbital)
!! Nbeta_elect=Number of beta electrons (N/2 for spin compensated systems)
!! Nalpha_elect=Number of beta electrons (N/2 for spin compensated systems)
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine rdm_init(RDMd,HighSpin,MSpin,INOF,Ista,NBF_occ,Nfrozen,Npairs,&
&  Ncoupled,Nbeta_elect,Nalpha_elect)
!Arguments ------------------------------------
!scalars
 logical,intent(in)::HighSpin
 integer,intent(in)::INOF,MSpin,Ista
 integer,intent(in)::NBF_occ,Nfrozen,Npairs,Ncoupled
 integer,intent(in)::Nbeta_elect,Nalpha_elect
 type(rdm_t),intent(inout)::RDMd
!Local variables ------------------------------
!scalars
!arrays

!************************************************************************

 RDMd%HighSpin=HighSpin
 RDMd%INOF=INOF
 RDMd%MSpin=MSpin
 RDMd%Ista=Ista
 RDMd%Nfrozen=Nfrozen
 RDMd%Nbeta_elect=Nbeta_elect
 RDMd%Nalpha_elect=Nalpha_elect
 RDMd%NBF_occ=NBF_occ
 RDMd%Ncoupled=Ncoupled
 RDMd%Npairs=Npairs
 RDMd%Nsingleocc=Nalpha_elect-Nbeta_elect
 RDMd%Npairs_p_sing=RDMd%Npairs+RDMd%Nsingleocc 
 RDMd%NBF_ldiag=RDMd%NBF_occ*(RDMd%NBF_occ+1)/2
 RDMd%Ngammas=RDMd%Ncoupled*RDMd%Npairs
 allocate(RDMd%DM2_J(RDMd%NBF_occ*RDMd%NBF_occ),RDMd%DM2_K(RDMd%NBF_occ*RDMd%NBF_occ)) 
 allocate(RDMd%Docc_gamma(RDMd%NBF_occ*RDMd%Ngammas)) 
 allocate(RDMd%DDM2_gamma_J(RDMd%NBF_occ*RDMd%NBF_occ*RDMd%Ngammas))
 allocate(RDMd%DDM2_gamma_K(RDMd%NBF_occ*RDMd%NBF_occ*RDMd%Ngammas)) 
 allocate(RDMd%OCC(RDMd%NBF_occ))
 allocate(RDMd%GAMMAs(RDMd%Ngammas))

end subroutine rdm_init
!!***

!!***
!!****f* DoNOF/rdm_free
!! NAME
!! rdm_free
!!
!! FUNCTION
!!  Free allocated arrays of the data type rdm_t 
!!
!! INPUTS
!!  
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine rdm_free(RDMd)
!Arguments ------------------------------------
!scalars
 class(rdm_t),intent(inout)::RDMd
!Local variables ------------------------------
!scalars
!arrays

!************************************************************************

 deallocate(RDMd%GAMMAs)
 deallocate(RDMd%OCC)
 deallocate(RDMd%DM2_J,RDMd%DM2_K) 
 deallocate(RDMd%Docc_gamma) 
 deallocate(RDMd%DDM2_gamma_J)
 deallocate(RDMd%DDM2_gamma_K) 

end subroutine rdm_free
!!***



end module m_rdmd
!!***

