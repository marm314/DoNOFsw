!!****m* DoNOF/m_optocc
!! NAME
!!  m_optocc
!!
!! FUNCTION
!!  Module prepared to compute occ optimization for a fixed hCOREpp and ERIs
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
 use m_cg
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
!!  imethocc=Method to use for the occ. optimization (default = 0, i.e. conjugate gradients)
!!  hCOREpp=diagonal part of the One-body integrals (h_pp) 
!!  ERI_J=Lower triangular part of the J_pq matrix
!!  ERI_K=Lower triangular part of the K_pq matrix
!!  Vnn=Nuclear-nuclear rep. energy
!!
!! OUTPUT
!!  Energy=Sum of nuclear and electronic energy
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine opt_occ(iter,imethod,RDMd,Vnn,Energy,hCOREpp,ERI_J,ERI_K) 
!Arguments ------------------------------------
!scalars
 integer,intent(inout)::iter
 integer,intent(in)::imethod
 double precision,intent(in)::Vnn
 double precision,intent(inout)::Energy
 type(rdm_t),intent(inout)::RDMd
!arrays
 double precision,dimension(RDMd%NBF_occ),intent(in)::hCOREpp
 double precision,dimension(RDMd%NBF_ldiag),intent(in)::ERI_J,ERI_K 
!Local variables ------------------------------
!scalars
 logical::diagco
 integer,parameter::msave=7,nextv=47,nfcall=6,nfgcal=7,g=28,toobig=2,vneed=4
 integer::iflag,ig,icall,icall1,Mtosave,Nwork,Nwork2
 double precision::eps,xtol
!arrays
 integer,dimension(2)::info_print
 integer,dimension(:),allocatable::iWork
 double precision,dimension(:),allocatable::GAMMAs,Grad_GAMMAs,diag,Work,Work2
!************************************************************************

 Energy=0.0d0
 allocate(GAMMAs(RDMd%Ngammas),GRAD_GAMMAs(RDMd%Ngammas))
 GRAD_GAMMAs=0.0d0
 if(iter==0) then 
  GAMMAs=0.785398163       ! Perturbed occ. numbers (i.e pi/4) -> occ(i<Fermi level) = 0.75
 else
  GAMMAs=RDMd%GAMMAs_old   ! Read from previous run
 endif

 icall=0
 if(imethod==1) then ! Conjugate gradients
  write(*,'(a)') 'Calling CG to optimize occ. numbers'
  Nwork=60; Nwork2=71+RDMd%Ngammas*(RDMd%Ngammas+15)/2; 
  allocate(iWork(Nwork),Work(RDMd%Ngammas),Work2(Nwork2))
  iWork=0; Work=0.1d0;   
  if (iWork(1)==0) call deflt(2,iWork, Nwork, Nwork2, Work2)
  iflag = iWork(1)
  if (iflag == 12 .or. iflag == 13) iWork(vneed) = iWork(vneed) + RDMd%Ngammas
  if (iflag == 14) goto 10
  if (iflag > 2 .and. iflag < 12) goto 10
  ig = 1
  if (iflag == 12) iWork(1) = 13
  goto 20

10   ig = iWork(g)

20   call sumit(Work, Energy, Work2(ig), iWork, Nwork, Nwork2, RDMd%Ngammas, Work2, GAMMAs)
  if(iWork(1)-2<0) then 
   goto 30
  elseif(iWork(1)-2==0) then
   goto 40
  else
   goto 50
  endif

30   icall1 = iWork(nfcall)
  call calc_E_occ(RDMd,GAMMAs,Energy,hCOREpp,ERI_J,ERI_K)
  icall=icall+1
  if (icall1 <= 0) iWork(toobig) = 1
  goto 20

40   call calc_Grad_occ(RDMd,Grad_GAMMAs,hCOREpp,ERI_J,ERI_K)
  Work2(ig:ig+RDMd%Ngammas)=Grad_GAMMAs(1:RDMd%Ngammas)
  goto 20

50   if(iWork(1) /= 14) then
        goto 60
     end if
!
!  Storage allocation
!
  iWork(g) = iWork(nextv)
  iWork(nextv) = iWork(g) + RDMd%Ngammas
  if(iflag /= 13) goto 10

60 deallocate(iWork,Work,Work2)

 else ! LBFGS
  write(*,'(a)') 'Calling LBFGS to optimize occ. numbers'
  Nwork=RDMd%Ngammas*(2*msave+1)+2*msave
  Mtosave=5; info_print(1)= -1; info_print(2)= 0; diagco= .false.;
  eps= 1.0d-5; xtol= 1.0d-16; icall=0; iflag=0;
  allocate(Work(Nwork),diag(RDMd%Ngammas))
  do
   call calc_E_occ(RDMd,GAMMAs,Energy,hCOREpp,ERI_J,ERI_K)
   call calc_Grad_occ(RDMd,Grad_GAMMAs,hCOREpp,ERI_J,ERI_K)
   call LBFGS(RDMd%Ngammas,Mtosave,GAMMAs,Energy,GRAD_GAMMAs,diagco,diag,info_print,eps,xtol,Work,iflag)
   if(iflag.le.0) exit
    icall=icall+1
!  We allow at most 2000 evaluations of Energy and Gradient
    if(icall.gt.2000) exit
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --       
  enddo
  deallocate(Work,diag)
 endif
 
 iter=iter+1
 if(iter>0) RDMd%GAMMAs_old=GAMMAs

 call calc_E_occ(RDMd,GAMMAs,Energy,hCOREpp,ERI_J,ERI_K)
 write(*,'(a,f15.6,a,i6,a)') 'Occ. optimized energy= ',Energy+Vnn,' after ',icall,' iter.'
 write(*,'(a)') ' '
 
 if(icall.gt.2000) write(*,'(a)') 'Warning! Max. number of iterations reached in occ. optimization'

 deallocate(GAMMAs,Grad_GAMMAs)

end subroutine opt_occ
!!***

end module m_optocc 
!!***

