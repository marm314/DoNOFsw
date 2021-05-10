!!****m* DoNOF/m_optorb
!! NAME
!!  m_optorb
!!
!! FUNCTION
!!  Module prepared to compute oorb optimization for a fixed OCCs and 2-RDM
!!
!! PARENTS
!!  m_noft_driver
!!
!! CHILDREN
!!  m_elag
!!  m_diagF
!!  m_e_grad_occ
!!
!! SOURCE
module m_optorb

 use m_rdmd
 use m_integd
 use m_elag
 use m_diagF
 use m_E_grad_occ

 implicit none

 private :: lambda_conv
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
!!  iter=Number of global iteration
!!  Vnn=Nuclear-nuclear rep. energy
!!
!! OUTPUT
!!  Energy=Sum of nuclear and electronic energy
!!  hCORE=One-body integrals (h_pq) in the optimized MO basis 
!!  ERI_J=Lower triangular part of the J_pq matrix in the optimized MO basis
!!  ERI_K=Lower triangular part of the K_pq matrix in the optimized MO basis
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine opt_orb(iter,imethod,ELAGd,RDMd,INTEGd,tol_dif_Lambda,Vnn,Energy,NO_COEF,mo_ints) 
!Arguments ------------------------------------
!scalars
 integer,intent(in)::iter,imethod
 double precision,intent(in)::Vnn,tol_dif_Lambda
 double precision,intent(inout)::Energy
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
!arrays
 double precision,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::NO_COEF
!Local variables ------------------------------
!scalars
 logical::convLambda
 integer::icall
!arrays
!************************************************************************

 Energy=0.0d0; convLambda=.false.;
 
 icall=0
 do
  icall=icall+1
  call mo_ints(NO_COEF,INTEGd%hCORE,INTEGd%ERImol)
  call ELAGd%build(RDMd,INTEGd,RDMd%DM2_J,RDMd%DM2_K)
  convLambda=lambda_conv(ELAGd,RDMd,tol_dif_Lambda)
  if(convLambda) exit  
  if(imethod==1) then ! Build F matrix for iterative diagonalization
   call diagF_to_coef(icall,ELAGd,RDMd,NO_COEF)
  else                ! Use Newton method to compute new COEFs
   
  endif
! We allow at most 2000 evaluations of Energy and Gradient
  if(icall.gt.2000) exit ! MAU
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
 enddo
 

 call INTEGd%eritoeriJK(RDMd%NBF_occ)
 call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K)
 write(*,'(a,f15.6,a,i6,a)') 'Orb. optimized energy= ',Energy+Vnn,' after ',icall,' iter.'
 
 if(icall.gt.2000) write(*,'(a)') 'Warning! Max. number of iterations reached in orb. optimization'

end subroutine opt_orb
!!***

!!***
!!****f* DoNOF/lambda_conv
!! NAME
!!  lambda_conv
!!
!! FUNCTION
!!  Check if the Lambda matrix already fulfils the condition Lambda_pq-Lambda_qp^* <= tol_dif_lambda.
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

function lambda_conv(ELAGd,RDMd,tol_dif_Lambda) result(converg_lamb)
!Arguments ------------------------------------
!scalars
 logical::converg_lamb
 double precision,intent(in)::tol_dif_Lambda
 type(elag_t),intent(in)::ELAGd
 type(rdm_t),intent(in)::RDMd
!arrays
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1
!arrays
!************************************************************************

 converg_lamb=.true.
 
 do iorb=1,RDMd%NBF_tot
  do iorb1=1,iorb-1
   if(dabs( ELAGd%Lambdas(iorb,iorb1)-ELAGd%Lambdas(iorb1,iorb) )>=tol_dif_Lambda .and. converg_lamb) then
    converg_lamb=.false.
   endif
  enddo
 enddo

end function lambda_conv
!!***

end module m_optorb
!!***
