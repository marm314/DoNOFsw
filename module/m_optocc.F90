!!****m* DoNOF/m_optocc
!! NAME
!!  m_optocc
!!
!! FUNCTION
!!  Module prepared to compute occ optimization for a fixed hCORE and ERIs
!!
!! PARENTS
!!  m_noft_driver
!!
!! CHILDREN
!!  m_E_grad_occ
!!
!! SOURCE
module m_optocc

 use m_rdmd
 use m_lbfgs
 use m_E_grad_occ

 implicit none

!!private :: 
!!***

 public :: opt_occ 
!!***

contains

!!***
!!****f* DoNOF/opt_occ
!! NAME
!!  opt_occ
!!
!! FUNCTION
!!  Call the CG of LBFGS subroutines for occ optimization 
!!
!! INPUTS
!!  hCORE=DiagoRDMd%Nalpha_electl part of the One-body integrals (h_pp) 
!!  ERI_J=Lower triangular part of the J_pq matrix
!!  ERI_K=Lower triangular part of the K_pq matrix
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine opt_occ(RDMd,hCORE,ERI_J,ERI_K) 
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
 double precision,dimension(RDMd%NBF_occ),intent(in)::hCORE
 double precision,dimension(RDMd%NBF_ldiag),intent(in)::ERI_J,ERI_K 
!Local variables ------------------------------
!scalars
 logical::diagco
 integer,parameter::msave=7
 integer::iflag,icall,Mtosave,Nwork
 double precision::Energy,eps,xtol
!arrays
 integer,dimension(2)::info_print
 double precision,dimension(:),allocatable::GAMMAs,Grad_GAMMAs,DIAG,Work
!************************************************************************

 Nwork=RDMd%Ngammas*(2*MSAVE+1)+2*msave
 allocate(GAMMAs(RDMd%Ngammas),GRAD_GAMMAs(RDMd%Ngammas))
 GAMMAs=0.1d0; GRAD_GAMMAs=0.0d0;


 Mtosave=5; info_print(1)= -1; info_print(2)= 0; diagco= .false.;
 eps= 1.0d-5; xtol= 1.0d-16; icall=0; iflag=0;
 allocate(Work(Nwork),DIAG(RDMd%Ngammas))
 do
  call calc_E_occ(RDMd,GAMMAs,Energy,hCORE,ERI_J,ERI_K)
  call calc_Grad_occ(RDMd,Grad_GAMMAs,hCORE,ERI_J,ERI_K)
  call LBFGS(RDMd%Ngammas,Mtosave,GAMMAs,Energy,GRAD_GAMMAs,diagco,DIAG,info_print,eps,xtol,Work,iflag)
  if(iflag.le.0) exit
   icall=icall + 1
!  We allow at most 2000 evaluations of Energy and Gradient
   if(icall.gt.2000) exit
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
 enddo
 deallocate(Work,DIAG)



 call calc_E_occ(RDMd,GAMMAs,Energy,hCORE,ERI_J,ERI_K)
 write(*,'(a,f15.6)') 'Optimized energy= ',Energy

 deallocate(GAMMAs,Grad_GAMMAs)

end subroutine opt_occ
!!***

end module m_optocc 
!!***

