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
 use m_e_grad_occ

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
!!  imethod=Method to use 1 -> F diag. 
!!  Vnn=Nuclear-nuclear rep. energy
!!  tol_dif_Lambda=Tolerance used to consider that F_pq is diagonal in the NO_COEF basis
!!  mo_ints=External subroutine that computes the hCORE and ERI integrals
!!
!! OUTPUT
!!  Energy=Sum of nuclear and electronic energy
!!  NO_COEF=Nat. orb. coefficients
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
 double precision::sumdiff,maxdiff
!arrays
!************************************************************************

 Energy=0.0d0; convLambda=.false.;
 
 icall=0
 call mo_ints(NO_COEF,INTEGd%hCORE,INTEGd%ERImol)
 do
  call ELAGd%build(RDMd,INTEGd,RDMd%DM2_J,RDMd%DM2_K)
  if(imethod==1) then ! Build F matrix for iterative diagonalization
   call lambda_conv(ELAGd,RDMd,tol_dif_Lambda,convLambda,sumdiff,maxdiff)
   if(convLambda) exit ! The NO_COEF and RDMs are already the solution =)
   call diagF_to_coef(iter,icall,maxdiff,ELAGd,RDMd,NO_COEF) ! Build new NO_COEF and set icall=icall+1
   if((iter==0).and.ELAGd%diagLpL_done) exit  ! We did Diag[(Lambda+Lambda)/2]. -> Do only one icall iteration before the occ. opt.
  else                ! Use Newton method to compute new COEFs
   
  endif
! Build all integrals in the new NO_COEF basis
  call mo_ints(NO_COEF,INTEGd%hCORE,INTEGd%ERImol) 
! We allow at most 50 evaluations of Energy and Gradient
  if(icall>50) exit 
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
 enddo
 
 ! Build ERI_J and ERI_K in NO_COEF basis before the occ. optimization 
 call INTEGd%eritoeriJK(RDMd%NBF_occ)
 call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K)
 write(*,'(a,f15.6,a,i6,a)') 'Orb. optimized energy= ',Energy+Vnn,' after ',icall,' iter.'
 
 if(icall>50) write(*,'(a)') 'Warning! Max. number of iterations (50) reached in orb. optimization'

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

subroutine lambda_conv(ELAGd,RDMd,tol_dif_Lambda,converg_lamb,sumdiff,maxdiff)
!Arguments ------------------------------------
!scalars
 logical,intent(inout)::converg_lamb
 double precision,intent(in)::tol_dif_Lambda
 double precision,intent(inout)::sumdiff,maxdiff
 type(elag_t),intent(in)::ELAGd
 type(rdm_t),intent(in)::RDMd
!arrays
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1
 double precision::diff
!arrays
!************************************************************************

 converg_lamb=.true.; sumdiff=0.0d0; maxdiff=0.0d0;
 
 do iorb=1,RDMd%NBF_tot
  do iorb1=1,RDMd%NBF_tot
   diff=dabs( ELAGd%Lambdas(iorb,iorb1)-ELAGd%Lambdas(iorb1,iorb) )
   sumdiff=sumdiff+diff
   if((diff>=tol_dif_Lambda) .and. converg_lamb) then
    converg_lamb=.false.
   endif
   if(diff>maxdiff) then
    maxdiff=diff
   endif
  enddo
 enddo

end subroutine lambda_conv
!!***

end module m_optorb
!!***
