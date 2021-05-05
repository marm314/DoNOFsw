!!****m* DoNOF/m_elag
!! NAME
!!  m_elag
!!
!! FUNCTION
!!  Module to build the Lagrange multipliers Lambda 
!!
!! PARENTS
!!  m_optorb
!!
!! CHILDREN
!!  m_rdmd
!!
!! SOURCE

module m_elag

 use m_rdmd
 use m_integd

 implicit none

!!private :: 
!!***

 public :: build_elag
!!***

contains
!!***

!!****f* DoNOF/build_elag
!! NAME
!! build_elag
!!
!! FUNCTION
!!  Build the Lagrange multipliers Lambda matrix. Nothe that the electron rep. integrals are given in DoNOF format
!!
!! INPUTS
!!  INTEGd=Object containg all integrals
!!  RDMd=Object containg all required variables whose arrays are properly updated
!!
!! OUTPUT
!!  Lambdas=Matrix build with the Lagrange multipliers Lambda_pq
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine build_elag(RDMd,INTEGd,Lambdas,DM2_J,DM2_K)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(in)::RDMd
 type(integ_t),intent(in)::INTEGd
!arrays
 double precision,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::Lambdas
 double precision,dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1
!arrays
!************************************************************************

 Lambdas=0.0d0

 do iorb=1,RDMd%NBF_occ
  Lambdas(iorb,:)=RDMd%OCC(iorb)*INTEGd%hCORE(:,iorb)                                  ! Init: Lambda_pq = n_p hCORE_qp
  Lambdas(iorb,:)=Lambdas(iorb,:)+RDMd%OCC(iorb)*INTEGd%ERImol(:,iorb,iorb,iorb)       ! any<->iorb,iorb<->iorb
  do iorb1=1,RDMd%NBF_occ
   Lambdas(iorb,:)=Lambdas(iorb,:)+DM2_J(iorb,iorb1)*INTEGd%ERImol(:,iorb1,iorb1,iorb) ! any<->iorb,iorb1<->iorb1
   Lambdas(iorb,:)=Lambdas(iorb,:)-DM2_K(iorb,iorb1)*INTEGd%ERImol(:,iorb1,iorb,iorb1) ! any<->iorb1,iorb1<->iorb
  enddo
 enddo 

 Lambdas=2.0d0*Lambdas

end subroutine build_elag
!!***

end module m_elag
!!***
