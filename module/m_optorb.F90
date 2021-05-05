!!****m* DoNOF/m_optorb
!! NAME
!!  m_optorb
!!
!! FUNCTION
!!  Module prepared to compute occ optimization for a fixed hCOREpp and ERIs
!!
!! PARENTS
!!  m_noft_driver
!!
!! CHILDREN
!!
!! SOURCE
module m_optorb

 use m_rdmd
 use m_integd
 use m_elag
 use m_E_grad_occ

 implicit none

!!private :: 
!!***

 public :: opt_orb 
!!***

contains

!!***
!!****f* DoNOF/opt_orb
!! NAME
!!  opt_orb
!!
!! FUNCTION
!!  Call the F-matrix method or Newton (Hessian) for orb optimization 
!!
!! INPUTS
!!  Vnn=Nuclear-nuclear rep. energy
!!
!! OUTPUT
!!  Energy=Sum of nuclear and electronic energy
!!  hCOREpp=diagonal part of the One-body integrals (h_pp) in the optimized MO basis 
!!  ERI_J=Lower triangular part of the J_pq matrix in the optimized MO basis
!!  ERI_K=Lower triangular part of the K_pq matrix in the optimized MO basis
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine opt_orb(RDMd,INTEGd,Vnn,Energy,MO_COEF,mo_ints) 
!Arguments ------------------------------------
!scalars
 double precision,intent(in)::Vnn
 double precision,intent(inout)::Energy
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
!arrays
 double precision,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::MO_COEF
!Local variables ------------------------------
!scalars
 integer::icall
!arrays
 double precision,dimension(:,:),allocatable::Lambdas
!************************************************************************

 allocate(Lambdas(RDMd%NBF_tot,RDMd%NBF_tot))
 Energy=0.0d0
 
 icall=0
 do
  icall=icall+1
  call mo_ints(MO_COEF,INTEGd%hCORE,INTEGd%ERImol)
  call build_elag(RDMd,INTEGd,Lambdas,RDMd%DM2_J,RDMd%DM2_K)

! We allow at most 2000 evaluations of Energy and Gradient
  if(icall.gt.0) exit ! MAU
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
 enddo
 
 deallocate(Lambdas)

 call INTEGd%htohpp(RDMd%NBF_occ)
 call INTEGd%eritoeriJK(RDMd%NBF_occ)
 call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy,INTEGd%hCOREpp,INTEGd%ERI_J,INTEGd%ERI_K)
 write(*,'(a,f15.6,a,i6,a)') 'Orb. optimized energy= ',Energy+Vnn,' after ',icall,' iter.'
 
 if(icall.gt.2000) write(*,'(a)') 'Warning! Max. number of iterations reached in orb. optimization'

end subroutine opt_orb
!!***

end module m_optorb
!!***
