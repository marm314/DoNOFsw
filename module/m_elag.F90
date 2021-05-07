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

 public :: build_elag,diag_lambda_ekt
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
!!  RDMd=Object containg all required variables whose arrays are properly updated
!!  INTEGd=Object containg all integrals
!!
!! OUTPUT
!!  RDMd%Lambdas=Matrix build with the Lagrange multipliers Lambda_pq
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine build_elag(RDMd,INTEGd,DM2_J,DM2_K)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(in)::INTEGd
!arrays
 double precision,dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1
!arrays
!************************************************************************

 RDMd%Lambdas=0.0d0

 do iorb=1,RDMd%NBF_occ
  RDMd%Lambdas(iorb,:)=RDMd%OCC(iorb)*INTEGd%hCORE(:,iorb)                                        ! Init: Lambda_pq = n_p hCORE_qp
  RDMd%Lambdas(iorb,:)=RDMd%Lambdas(iorb,:)+RDMd%OCC(iorb)*INTEGd%ERImol(:,iorb,iorb,iorb)        ! any<->iorb,iorb<->iorb
  do iorb1=1,RDMd%NBF_occ
   if(iorb/=iorb1) then
    RDMd%Lambdas(iorb,:)=RDMd%Lambdas(iorb,:)+DM2_J(iorb,iorb1)*INTEGd%ERImol(:,iorb1,iorb1,iorb) ! any<->iorb,iorb1<->iorb1
    RDMd%Lambdas(iorb,:)=RDMd%Lambdas(iorb,:)-DM2_K(iorb,iorb1)*INTEGd%ERImol(:,iorb1,iorb,iorb1) ! any<->iorb1,iorb1<->iorb
   endif
  enddo
 enddo 

 ! TODO 
 if(RDMd%Nsingleocc>0) write(*,'(a)') 'Error! The Lambda_pq matrix construction is not implemented for Nsingleocc>0'

 !RDMd%Lambdas=2.0d0*RDMd%Lambdas ! We only need half for 'alpha' orbs to define gradients

end subroutine build_elag
!!***

!!****f* DoNOF/diag_lambda_ekt
!! NAME
!! diag_lambda_ekt
!!
!! FUNCTION
!!  Diagonalize the Lagrange multipliers Lambda matrix (produce either the 'canonical orbitals' or EKT). 
!!
!! INPUTS
!!  RDMd=Object containg all required variables whose arrays are properly updated
!!  INTEGd=Object containg all integrals
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine diag_lambda_ekt(RDMd,NO_COEF,ekt)
!Arguments ------------------------------------
!scalars
 logical,optional::ekt
 type(rdm_t),intent(in)::RDMd
!arrays
 double precision,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::NO_COEF
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,lwork,info
 double precision::sqrt_occ_iorb,sqrt_occ_iorb1,tol6=1d-6
!arrays
 character(len=10)::coef_file
 double precision,dimension(:),allocatable::Eigval,Eigval_occ,Work
 double precision,dimension(:,:),allocatable::Eigvec,CANON_COEF
!************************************************************************

 allocate(Eigvec(RDMd%NBF_tot,RDMd%NBF_tot),Eigval(RDMd%NBF_tot),Work(1))
 allocate(Eigval_occ(RDMd%NBF_occ))
 
 Eigvec=RDMd%Lambdas

 if(present(ekt)) then
  do iorb=1,RDMd%NBF_tot
   do iorb1=1,RDMd%NBF_tot
    if(iorb<=RDMd%NBF_occ.and.iorb1<=RDMd%NBF_occ) then
     sqrt_occ_iorb =DSQRT(RDMd%OCC(iorb))
     sqrt_occ_iorb1=DSQRT(RDMd%OCC(iorb1))
     if((DABS(sqrt_occ_iorb)>tol6).and.(DABS(sqrt_occ_iorb1)>tol6)) then
      Eigvec(iorb,iorb1)=Eigvec(iorb,iorb1)/(sqrt_occ_iorb*sqrt_occ_iorb1)
     else
      Eigvec(iorb,iorb1)=0.0d0
     endif
    else
     Eigvec(iorb,iorb1)=0.0d0
    endif
   enddo
  enddo
 endif

 lwork=-1 
 call DSYEV('V','L',RDMd%NBF_tot,Eigvec,RDMd%NBF_tot,Eigval,Work,lwork,info)
 lwork= nint(Work(1))

 if(info==0) then
  deallocate(Work)
  allocate(Work(lwork)) 
  call DSYEV('V','L',RDMd%NBF_tot,Eigvec,RDMd%NBF_tot,Eigval,Work,lwork,info)
 endif

 ! Print final eigenvalues
 write(*,'(a)') ' '
 if(present(ekt)) then
  Eigval=-Eigval
  write(*,'(a)') 'EKT ionization potentials'
 else
  coef_file='CANON_COEF'
  allocate(CANON_COEF(RDMd%NBF_tot,RDMd%NBF_tot))
  CANON_COEF=matmul(NO_COEF,Eigvec)
  call RDMd%print_orbs(CANON_COEF,coef_file)
  deallocate(CANON_COEF)
  write(*,'(a)') 'Canonical orbital eigenvalues'
 endif

 Eigval_occ(1:RDMd%NBF_occ)=Eigval(1:RDMd%NBF_occ)
 do iorb=1,(RDMd%NBF_occ/10)*10,10
  write(*,'(f12.6,9f11.6)') Eigval_occ(iorb:iorb+9)
 enddo
 iorb=(RDMd%NBF_occ/10)*10+1
 write(*,'(f12.6,*(f11.6))') Eigval_occ(iorb:)
 write(*,'(a)') ' '

  
 deallocate(Eigvec,Work,Eigval,Eigval_occ)

end subroutine diag_lambda_ekt
!!***

end module m_elag
!!***
