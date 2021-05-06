!!****m* DoNOF/m_integd
!! NAME
!! basic integrals variables
!!
!! FUNCTION
!! This module contains definitions of the arrays containing all integrals used by NOFT module.
!!
!! COPYRIGHT
!! This file is distributed under the terms of the
!! GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! SOURCE

module m_integd

 implicit none

!!****t* m_integd/integ_t
!! NAME
!! integ_t
!!
!! FUNCTION
!! Datatype storing integral arrays needed
!!
!! SOURCE

 type,public :: integ_t

! arrays 
 double precision,dimension(:),allocatable::hCOREpp,ERI_J,ERI_K
 double precision,dimension(:,:),allocatable::hCORE,Lambdas
 double precision,dimension(:,:,:,:),allocatable::ERImol

 contains 
   procedure :: free => integ_free
   ! Destructor.

   procedure :: htohpp => hcore_to_hcorepp
   ! hCORE to hCOREpp.

   procedure :: eritoeriJK => eri_to_eriJK  
   ! ERImol to ERI_J and ERI_K.

 end type integ_t

 public :: integ_init     ! Main creation method.
!!***

CONTAINS  !==============================================================================

!!***
!!****f* DoNOF/integ_init
!! NAME
!! integ_init
!!
!! FUNCTION
!!  Initialize the data type integ_t 
!!
!! INPUTS
!! NBF_tot=Number of total orbitals
!! NBF_occ=Number of orbitals that are occupied
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine integ_init(Integd,NBF_tot,NBF_occ)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::NBF_tot,NBF_occ
 type(integ_t),intent(inout)::Integd
!Local variables ------------------------------
!scalars
 integer::NBF_ldiag
!arrays
!************************************************************************

 NBF_ldiag=NBF_occ*(NBF_occ+1)/2
 allocate(Integd%hCOREpp(NBF_occ),Integd%ERI_J(NBF_ldiag),Integd%ERI_K(NBF_ldiag))
 allocate(Integd%hCORE(NBF_tot,NBF_tot),Integd%Lambdas(NBF_tot,NBF_tot))
 allocate(Integd%ERImol(NBF_tot,NBF_tot,NBF_tot,NBF_tot))

end subroutine integ_init
!!***

!!***
!!****f* DoNOF/integ_free
!! NAME
!! integ_free
!!
!! FUNCTION
!!  Free allocated arrays of the data type integ_t 
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

subroutine integ_free(Integd)
!Arguments ------------------------------------
!scalars
 class(integ_t),intent(inout)::Integd
!Local variables ------------------------------
!scalars
!arrays
!************************************************************************

 deallocate(Integd%Lambdas) 
 deallocate(Integd%hCORE) 
 deallocate(Integd%hCOREpp) 
 deallocate(Integd%ERImol) 
 deallocate(Integd%ERI_J) 
 deallocate(Integd%ERI_K) 

end subroutine integ_free
!!***

!!***
!!****f* DoNOF/hcore_to_hcorepp
!! NAME
!! hcore_to_hcorepp
!!
!! FUNCTION
!!  Get hCOREpp from hCORE 
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

subroutine hcore_to_hcorepp(Integd,NBF_occ)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::NBF_occ
 class(integ_t),intent(inout)::Integd
!Local variables ------------------------------
!scalars
 integer::iorb
!arrays
!************************************************************************

 do iorb=1,NBF_occ
  Integd%hCOREpp(iorb)=Integd%hCORE(iorb,iorb)
 enddo

end subroutine hcore_to_hcorepp
!!***

!!***
!!****f* DoNOF/eri_to_eriJK
!! NAME
!! eri_to_eriJK
!!
!! FUNCTION
!!  Get ERI_J, ERI_K from ERI
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

subroutine eri_to_eriJK(Integd,NBF_occ)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::NBF_occ
 class(integ_t),intent(inout)::Integd
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2
!arrays
!************************************************************************

 iorb2=1
 do iorb=1,NBF_occ
  do iorb1=1,iorb
   Integd%ERI_J(iorb2)=Integd%ERImol(iorb,iorb1,iorb1,iorb) ! J in DoNOF
   Integd%ERI_K(iorb2)=Integd%ERImol(iorb,iorb1,iorb,iorb1) ! K in DoNOF
   iorb2=iorb2+1
  enddo
 enddo

end subroutine eri_to_eriJK
!!***

end module m_integd
!!***
