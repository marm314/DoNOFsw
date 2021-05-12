!!****m* DoNOF/m_elag
!! NAME
!!  m_elag
!!
!! FUNCTION
!!  Module to build the Lagrange multipliers Lambda 
!!
!! COPYRIGHT
!! This file is distributed under the terms of the
!! GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
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
!!****t* m_rdmd/rdm_t
!! NAME
!! rdm_t
!!
!! FUNCTION
!! Datatype storing noft quantities and arrays needed
!!
!! SOURCE

 type,public :: elag_t

  logical::diagLpL=.true.        ! Do the diag. using (lambda+lambda)/2?
  logical::diagLpL_done=.false.  ! Did we use use (lambda+lambda)/2?
  integer::MaxScaling=0          ! Max scaling reductions employed to avoid divergence of diag[F]
  integer::itscale=1             ! Above this number of iterations we do MaxScaling=MaxScaling+1
  integer::itolLambda=4          ! Integer used to define 10**-itolLambda as threshold of Lambda_pq-Lambda_qp* convergence
  integer::itoldiis=3            ! Integer used to define 10**-itoldiis as threshold of DIIS trigger
  integer::idiis=0               ! Current DIIS iteration
  integer::ndiis=5               ! The number of iterations required to apply DIIS is ndiis+1
  integer::ndiis_array           ! Size of the arrays used in DIIS (ndiis+2)
  double precision::sumdiff_old  ! Old value of sum_pq |F_pq|  for p/=q 
! arrays 
  double precision,allocatable,dimension(:)::F_diag       ! F_pp (Diag. part of the F matrix)
  double precision,allocatable,dimension(:)::Coef_DIIS    ! DIIS coefs. used to build linear comb. of F matrices
  double precision,allocatable,dimension(:,:)::Lambdas    ! Lambda_pq (Lagrange multipliers matrix)
  double precision,allocatable,dimension(:,:)::DIIS_mat   ! DIIS matrix used to solve the system of eqs. DIIS_MAT*Coef_DIIS = (0 0 ... 0 1) 
  double precision,allocatable,dimension(:,:,:)::F_DIIS   ! F matrices used by DIIS

 contains 
   procedure :: free => elag_free
   ! Destructor.

   procedure :: build => build_elag
   ! Use integrals and the 1,2-RDM to build Lambdas matrix.

   procedure :: diag_lag => diag_lambda_ekt
   ! Diagonalize the matrix Lambdas (or divided by occ. numbers) to compute canonical orbs. or EKT.

   procedure :: clean_diis => wipeout_diis
   ! Set to ZERO all arrays employed by DIIS.

 end type elag_t

 public :: elag_init 
!!***

CONTAINS  !==============================================================================

!!***
!!****f* DoNOF/elag_init
!! NAME
!! elag_init
!!
!! FUNCTION
!!  Initialize the data type elag_t 
!!
!! INPUTS
!! NBF_tot=Number of total orbitals
!! diagLpL_in=Diagonalize 0.5 (Lambda+Lambda) for the first iteration?
!! itolLambda=Used as 10**-itolLambda to check for Lambda_pq-Lambda_qp* convergence
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine elag_init(ELAGd,NBF_tot,diagLpL_in,itolLambda_in,ndiis_in)
!Arguments ------------------------------------
!scalars
 logical,intent(in)::diagLpL_in
 integer,intent(in)::NBF_tot,itolLambda_in,ndiis_in
 type(elag_t),intent(inout)::ELAGd
!Local variables ------------------------------
!scalars
!arrays
!************************************************************************

 ELAGd%itolLambda=itolLambda_in
 ELAGd%diagLpL=diagLpL_in
 ELAGd%ndiis=ndiis_in
 ELAGd%ndiis_array=ELAGd%ndiis+2
 allocate(ELAGd%F_diag(NBF_tot))
 allocate(ELAGd%Lambdas(NBF_tot,NBF_tot)) 
 if(ELAGd%ndiis>0) then
  allocate(ELAGd%Coef_DIIS(ELAGd%ndiis_array))
  allocate(ELAGd%F_DIIS(ELAGd%ndiis_array,NBF_tot,NBF_tot))
  allocate(ELAGd%DIIS_mat(ELAGd%ndiis_array,ELAGd%ndiis_array)) 
 endif 
 
end subroutine elag_init
!!***

!!***
!!****f* DoNOF/elag_free
!! NAME
!! elag_free
!!
!! FUNCTION
!!  Free allocated arrays of the data type elag_t 
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

subroutine elag_free(ELAGd)
!Arguments ------------------------------------
!scalars
 class(elag_t),intent(inout)::ELAGd
!Local variables ------------------------------
!scalars
!arrays
!************************************************************************

 deallocate(ELAGd%F_diag) 
 deallocate(ELAGd%Lambdas) 
 if(ELAGd%ndiis>0) then
  deallocate(ELAGd%Coef_DIIS)
  deallocate(ELAGd%F_DIIS)
  deallocate(ELAGd%DIIS_mat) 
 endif 

end subroutine elag_free
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
!!  ELAGd%Lambdas=Matrix build with the Lagrange multipliers Lambda_pq
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine build_elag(ELAGd,RDMd,INTEGd,DM2_J,DM2_K)
!Arguments ------------------------------------
!scalars
 class(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(in)::INTEGd
!arrays
 double precision,dimension(RDMd%NBF_occ,RDMd%NBF_occ),intent(inout)::DM2_J,DM2_K
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1
 double precision::tol10=1.0d-10
!arrays
!************************************************************************

 ELAGd%Lambdas=0.0d0

 do iorb=1,RDMd%NBF_occ
  ELAGd%Lambdas(iorb,:)=RDMd%occ(iorb)*INTEGd%hCORE(:,iorb)                                         ! Init: Lambda_pq = n_p hCORE_qp
  ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+RDMd%occ(iorb)*INTEGd%ERImol(:,iorb,iorb,iorb)        ! any<->iorb,iorb<->iorb
  do iorb1=1,RDMd%NBF_occ
   if(iorb/=iorb1) then
    ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)+DM2_J(iorb,iorb1)*INTEGd%ERImol(:,iorb1,iorb1,iorb) ! any<->iorb,iorb1<->iorb1
    ELAGd%Lambdas(iorb,:)=ELAGd%Lambdas(iorb,:)-DM2_K(iorb,iorb1)*INTEGd%ERImol(:,iorb1,iorb,iorb1) ! any<->iorb1,iorb1<->iorb
   endif
  enddo
 enddo 
 !ELAGd%Lambdas=2.0d0*ELAGd%Lambdas ! We only need half for 'alpha' orbs to define gradients

 ! TODO 
 if(RDMd%Nsingleocc>0) write(*,'(a)') 'Error! The Lambda_pq matrix construction is not implemented for Nsingleocc>0'

 do iorb=1,RDMd%NBF_tot
  do iorb1=1,RDMd%NBF_tot
   if(dabs(ELAGd%Lambdas(iorb,iorb1))<tol10) ELAGd%Lambdas(iorb,iorb1)=0.0d0
  enddo
 enddo 

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
!!  ELAGd%Lambdas=Matrix containing the Lagrange multipliers Lambda_pq
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

subroutine diag_lambda_ekt(ELAGd,RDMd,NO_COEF,ekt)
!Arguments ------------------------------------
!scalars
 logical,optional::ekt
 class(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(in)::RDMd
!arrays
 double precision,dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(in)::NO_COEF
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,lwork,info
 double precision::sqrt_occ_iorb,sqrt_occ_iorb1,tol6=1d-6
!arrays
 character(len=10)::coef_file
 double precision,allocatable,dimension(:)::Eigval,Eigval_nocc,Work
 double precision,allocatable,dimension(:,:)::Eigvec,CANON_COEF
!************************************************************************

 allocate(Eigvec(RDMd%NBF_tot,RDMd%NBF_tot),Eigval(RDMd%NBF_tot),Work(1))
 allocate(Eigval_nocc(RDMd%NBF_occ))
 
 Eigvec=ELAGd%Lambdas

 if(present(ekt)) then
  do iorb=1,RDMd%NBF_tot
   do iorb1=1,RDMd%NBF_tot
    if(iorb<=RDMd%NBF_occ.and.iorb1<=RDMd%NBF_occ) then
     sqrt_occ_iorb =dsqrt(RDMd%occ(iorb))
     sqrt_occ_iorb1=dsqrt(RDMd%occ(iorb1))
     if((dabs(sqrt_occ_iorb)>tol6).and.(dabs(sqrt_occ_iorb1)>tol6)) then
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
 lwork=nint(Work(1))

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

 Eigval_nocc(1:RDMd%NBF_occ)=Eigval(1:RDMd%NBF_occ)
 do iorb=1,(RDMd%NBF_occ/10)*10,10
  write(*,'(f12.6,9f11.6)') Eigval_nocc(iorb:iorb+9)
 enddo
 iorb=(RDMd%NBF_occ/10)*10+1
 write(*,'(f12.6,*(f11.6))') Eigval_nocc(iorb:)
 write(*,'(a)') ' '

  
 deallocate(Eigvec,Work,Eigval,Eigval_nocc)

end subroutine diag_lambda_ekt
!!***

!!****f* DoNOF/wipeout_diis
!! NAME
!! wipeout_diis
!!
!! FUNCTION
!!  Build the Lagrange multipliers Lambda matrix. Nothe that the electron rep. integrals are given in DoNOF format
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

subroutine wipeout_diis(ELAGd)
!Arguments ------------------------------------
!scalars
 class(elag_t),intent(inout)::ELAGd
!arrays
!Local variables ------------------------------
!scalars
!arrays
!************************************************************************

 ELAGd%idiis=0
 if(ELAGd%ndiis>0) then
  ELAGd%Coef_DIIS=0.0d0
  ELAGd%F_DIIS=0.0d0
  ELAGd%DIIS_mat=0.0d0
 endif 

end subroutine wipeout_diis
!!***

end module m_elag
!!***
