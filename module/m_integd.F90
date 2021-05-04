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
 double precision,dimension(:,:),allocatable::hCORE
 double precision,dimension(:,:,:,:),allocatable::ERImol

 contains 
   procedure :: free => integ_free
   ! Destructor.

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
 allocate(Integd%hCORE(NBF_tot,NBF_tot))
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

 deallocate(Integd%hCORE) 
 deallocate(Integd%hCOREpp) 
 deallocate(Integd%ERImol) 
 deallocate(Integd%ERI_J) 
 deallocate(Integd%ERI_K) 

end subroutine integ_free

end module m_integd
!!***
