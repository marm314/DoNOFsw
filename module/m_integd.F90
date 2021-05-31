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
 double precision,allocatable,dimension(:)::ERI_J,ERI_K
 double precision,allocatable,dimension(:,:)::hCORE,Overlap
 double precision,allocatable,dimension(:,:,:,:)::ERImol

 contains 
   procedure :: free => integ_free
   ! Destructor.

   procedure :: eritoeriJK => eri_to_eriJK  
   ! ERImol to ERI_J and ERI_K.

   procedure :: print_int => print_ints 
   ! Print hCORE and ERImol integrals in their current status

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
!! Overlap_in=S_ao overlap matrix in atomic orbs.
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine integ_init(Integd,NBF_tot,NBF_occ,Overlap_in)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::NBF_tot,NBF_occ
 type(integ_t),intent(inout)::Integd
 double precision,dimension(NBF_tot,NBF_tot),intent(in)::Overlap_in
!Local variables ------------------------------
!scalars
 integer::NBF_ldiag
 double precision::totMEM
!arrays
!************************************************************************

 NBF_ldiag=NBF_occ*(NBF_occ+1)/2
 ! Calculate memory needed
 totMEM=2*NBF_ldiag+2*NBF_tot*NBF_tot+NBF_tot*NBF_tot*NBF_tot*NBF_tot
 totMEM=8*totMEM       ! Bytes
 totMEM=totMEM*1.0d-6  ! Bytes to Mb  
 if(totMEM>1.0d3) then     ! Mb to Gb
  write(*,'(a,f10.3,a)') 'Mem. required for storing INTEGd object ',totMEM*1.0d-3,' Gb'
 elseif(totMEM<1.0d0) then ! Mb to Kb
  write(*,'(a,f10.3,a)') 'Mem. required for storing INTEGd object ',totMEM*1.0d3,' Kb'
 else                      ! Mb
  write(*,'(a,f10.3,a)') 'Mem. required for storing INTEGd object ',totMEM,' Mb'
 endif
 ! Allocate arrays
 allocate(Integd%ERI_J(NBF_ldiag),Integd%ERI_K(NBF_ldiag))
 allocate(Integd%hCORE(NBF_tot,NBF_tot),Integd%Overlap(NBF_tot,NBF_tot))
 allocate(Integd%ERImol(NBF_tot,NBF_tot,NBF_tot,NBF_tot))
 Integd%Overlap=Overlap_in

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
 deallocate(Integd%Overlap) 
 deallocate(Integd%ERImol) 
 deallocate(Integd%ERI_J) 
 deallocate(Integd%ERI_K) 

end subroutine integ_free
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

!!***
!!****f* DoNOF/print_ints
!! NAME
!! print_ints
!!
!! FUNCTION
!!  Print hCORE and ERImol integrals in their current status 
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

subroutine print_ints(Integd,NBF_tot)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::NBF_tot
 class(integ_t),intent(in)::Integd
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iunit=312
 double precision::tol8=1.0d-8
!arrays
!************************************************************************
 
 ! Print ERImol
 open(unit=iunit,form='unformatted',file='ERImol')
 do iorb=1,NBF_tot
  do iorb1=1,NBF_tot
   do iorb2=1,NBF_tot
    do iorb3=1,NBF_tot
     if(dabs(INTEGd%ERImol(iorb,iorb1,iorb2,iorb3))>tol8) then
      write(iunit) iorb,iorb1,iorb,iorb1,INTEGd%ERImol(iorb,iorb1,iorb2,iorb3)
     endif
    enddo
   enddo
  enddo
 enddo 
 write(iunit) 0,0,0,0,0.0d0
 write(iunit) 0,0,0,0,0.0d0
 close(iunit)

 ! Print hCORE
 open(unit=iunit,form='unformatted',file='hCORE')
 do iorb=1,NBF_tot
  do iorb1=1,NBF_tot
   if(dabs(INTEGd%hCORE(iorb,iorb1))>tol8) then
    write(iunit) iorb,iorb1,INTEGd%hCORE(iorb,iorb1)
   endif
  enddo
 enddo
 write(iunit) 0,0,0.0d0
 close(iunit)

end subroutine print_ints
!!***

end module m_integd
!!***
