!!****m* DoNOF/m_diagf
!! NAME
!!  m_diagf
!!
!! FUNCTION
!!  Module prepared to compute oorb optimization for a fixed OCCs and 2-RDM
!!
!! PARENTS
!!  m_optorb
!!
!! CHILDREN
!!
!! SOURCE
module m_diagf

 use m_rdmd
 use m_elag

 implicit none

!!private :: 
!!***

 public :: diagF_to_coef
!!***

contains

!!***
!!****f* DoNOF/diagF_to_coef
!! NAME
!!  diagF_to_coef
!!
!! FUNCTION
!!  Build the F-matrix, diagonalize it and update the NO_COEF
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

subroutine diagF_to_coef(icall,ELAGd,RDMd,NO_COEF) 
!Arguments ------------------------------------
!scalars
 integer,intent(in)::icall
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(in)::RDMd
!arrays
 double precision,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::NO_COEF
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,lwork,info
!arrays
 double precision,dimension(:),allocatable::Work
 double precision,dimension(:,:),allocatable::Eigvec,New_NO_COEF
!************************************************************************
 
 allocate(New_NO_COEF(RDMd%NBF_tot,RDMd%NBF_tot),Eigvec(RDMd%NBF_tot,RDMd%NBF_tot),Work(1))

 if((icall==1).and.ELAGd%diagLpL) then
  ELAGd%diagLpL=.false. 
  do iorb=1,RDMd%NBF_tot 
   Eigvec(iorb,iorb)=ELAGd%Lambdas(iorb,iorb)
   do iorb1=1,iorb-1
    Eigvec(iorb,iorb1)=0.5d0*(ELAGd%Lambdas(iorb,iorb1)+ELAGd%Lambdas(iorb1,iorb))
    Eigvec(iorb1,iorb)=Eigvec(iorb,iorb1)
   enddo
  enddo  
 else
  do iorb=1,RDMd%NBF_tot 
   Eigvec(iorb,iorb)=ELAGd%Lambdas_diag(iorb)
   do iorb1=1,iorb-1
    Eigvec(iorb,iorb1)=ELAGd%Lambdas(iorb,iorb1)-ELAGd%Lambdas(iorb1,iorb)
    Eigvec(iorb1,iorb)=Eigvec(iorb,iorb1)
   enddo
  enddo  
  ! Scale or sth TODO
 endif

 lwork=-1
 call DSYEV('V','L',RDMd%NBF_tot,Eigvec,RDMd%NBF_tot,ELAGd%Lambdas_diag,Work,lwork,info)
 lwork=nint(Work(1))

 if(info==0) then
  deallocate(Work)
  allocate(Work(lwork))
  call DSYEV('V','L',RDMd%NBF_tot,Eigvec,RDMd%NBF_tot,ELAGd%Lambdas_diag,Work,lwork,info)
 endif

 New_NO_COEF=matmul(NO_COEF,Eigvec)
 NO_COEF(:,:)=New_NO_COEF(:,:)

 deallocate(New_NO_COEF,Eigvec,Work)

end subroutine diagF_to_coef
!!***

end module m_diagf
!!***
