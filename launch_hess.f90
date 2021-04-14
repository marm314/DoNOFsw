      SUBROUTINE LAUNCH_HESS(ELAG,COEF,RO,CJ12,CK12,AHCORE,ADIPx,ADIPy,ADIPz,IERI,ERI)
      USE PARCOM
      PARAMETER (ZERO=0.0d0)
      INTEGER::i,j,k,l,a,b,nbf2,nbf3,nbf4
      DOUBLE PRECISION::Epnof,Epnof_old
      INTEGER,DIMENSION(NIJKL)::IERI
      DOUBLE PRECISION,DIMENSION(NIJKL)::ERI
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AHCORE,ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ELAG,COEF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::OCC
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::HCOREmol
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::VEC
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:)::ERImol
!-----------------------------------------------------------------------
!     NCO:  Number of HF occupied MOs (OCC=1 in SD)
!     NVIR: Number of HF virtual  MOs (OCC=0 in SD)
!
!     NO1:  Number of inactive doubly occupied orbitals (OCC=1)
!     NDOC: Number of strongly occupied MOs
!     NCWO: Number of coupled weakly occupied MOs per strongly occupied
!     NCWO*NDOC: Active orbitals in the virtual subspace
!     NO0: Empty orbitals  (OCC=0)
!     NAC:  Dimension of the active natural orbital subspace
!
!     NB/NCO: Number of orbitals in the doubly-occupied space
!     NA: NB/NCO + Number of orbitals in the singly-occupied space (NSOC)
!
!           NCO + NSOC     |       NVIR          = NBF
!-----------------------------------------------------------------------
      nbf2=nbf*nbf
      nbf3=nbf*nbf2
      nbf4=nbf*nbf4
      allocate(OCC(nbf),VEC(nbf,nbf),HCOREmol(nbf,nbf),ERImol(nbf,nbf,nbf,nbf))
      HCOREmol=ZERO;ERImol=ZERO;
      OCC(:)=RO(:)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  These 5 lines prepare HCOREmol and ERImol   !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      VEC=COEF
      HCOREmol=matmul(transpose(VEC),matmul(AHCORE,VEC))
      CALL ERIC1c(ERImol,IERI,ERI,VEC,NBF)
      CALL ERIC23c(ERImol,VEC,NBF)
      CALL ERIC4c(ERImol,VEC,NBF)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) 'Hello world',IPNOF,THRESHEN,MAXIT
      do i=1,MAXIT
       
       Epnof=EN+calc_epnof(OCC,HCOREmol,ERImol,NBF,NDOC,NO1,NCWO,IPNOF) 
       write(*,*) epnof
       if(abs(Epnof-Epnof_old)>THRESHEN) goto 92

      enddo 
92    continue      
      deallocate(OCC,VEC,HCOREmol,ERImol)

      END SUBROUTINE LAUNCH_HESS
    
      function calc_epnof(OCC,HCOREmol,ERImol,nbf,nfl,ninact,ncoup,ipnof)
      implicit none
      integer::nbf,nfl,ninact,ncoup,ipnof
      double precision::calc_epnof
      double precision,dimension(nbf)::OCC
      double precision,dimension(nbf,nbf)::HCOREmol
      double precision,dimension(nbf,nbf,nbf,nbf)::ERImol
      integer::iorb1,iorb2,iorb3,forb,lorb
      double precision::horb1,horb2,horb3
      calc_epnof=0.0d0
      do iorb1=1,nfl
       if(iorb1>ninact) then
        horb1=1.0d0-OCC(iorb1)
        calc_epnof=calc_epnof+OCC(iorb1)*(2.0d0*HCOREmol(iorb1,iorb1)+ERImol(iorb1,iorb1,iorb1,iorb1)) 
        ! Intra [same for PNOF5, PNOF7, and PNOF7s]
        ! np (2 hpp + Jpp)
        forb=(nfl+1)+ncoup*(nfl-iorb1)
        lorb=(nfl+ncoup)+ncoup*(nfl-iorb1)
        do iorb2=forb,lorb
         calc_epnof=calc_epnof+OCC(iorb2)*(2.0d0*HCOREmol(iorb2,iorb2)+ERImol(iorb2,iorb2,iorb2,iorb2)) 
         ! - sqrt(np nq) Lpq [below Fl and above ]
         calc_epnof=calc_epnof-sqrt(OCC(iorb1)*OCC(iorb2))*(ERImol(iorb1,iorb2,iorb1,iorb2)+ERImol(iorb2,iorb1,iorb2,iorb1))
         ! + sqrt(np nq) Lpq [above Fl both]
         do iorb3=iorb2+1,lorb
          calc_epnof=calc_epnof+sqrt(OCC(iorb2)*OCC(iorb3))*(ERImol(iorb3,iorb2,iorb3,iorb2)+ERImol(iorb2,iorb3,iorb2,iorb3))
         enddo
         ! Inter PNOF5 ->  np nq (2 Jpq - Kpq)  [above Fl both]
         do iorb3=lorb+1,nbf
          calc_epnof=calc_epnof+(OCC(iorb2)*OCC(iorb3))&
                    & *(2.0d0*(ERImol(iorb2,iorb3,iorb3,iorb2)+ERImol(iorb3,iorb2,iorb2,iorb3))&
                    & -(ERImol(iorb2,iorb3,iorb2,iorb3)+ERImol(iorb3,iorb2,iorb3,iorb2)))
         enddo
         ! Inter PNOF7 ->  - sqrt(np hp nq hq) Lpq [Omegag and Omegaf, above Fl both]
         if(ipnof==7) then
          horb2=1.0d0-OCC(iorb2)
          do iorb3=lorb+1,nbf
           horb3=1.0d0-OCC(iorb3)
           calc_epnof=calc_epnof-sqrt(OCC(iorb2)*OCC(iorb3)*horb2*horb3)&
                     & *(ERImol(iorb2,iorb3,iorb2,iorb3)+ERImol(iorb3,iorb2,iorb3,iorb2))
          enddo
         endif
        enddo
        ! Inter PNOF5 -> np nq (2 Jpq - Kpq) [Omegag and Omegaf, below and below Fl]
        do iorb2=iorb1+1,nfl
         calc_epnof=calc_epnof+(OCC(iorb2)*OCC(iorb1))&
                   & *(2.0d0*(ERImol(iorb2,iorb1,iorb1,iorb2)+ERImol(iorb1,iorb2,iorb2,iorb1))&
                   & -(ERImol(iorb2,iorb1,iorb2,iorb1)+ERImol(iorb1,iorb2,iorb1,iorb2)))
        enddo
        ! Inter PNOF5 -> np nq (2 Jpq - Kpq) [Omegag and Omegaf, below and above Fl]
        do iorb2=nfl+1,forb-1
         calc_epnof=calc_epnof+(OCC(iorb2)*OCC(iorb1))&
                   & *(2.0d0*(ERImol(iorb2,iorb1,iorb1,iorb2)+ERImol(iorb1,iorb2,iorb2,iorb1))&
                   & -(ERImol(iorb2,iorb1,iorb2,iorb1)+ERImol(iorb1,iorb2,iorb1,iorb2)))
        enddo
        do iorb2=lorb+1,nbf
         calc_epnof=calc_epnof+(OCC(iorb2)*OCC(iorb1))&
                   & *(2.0d0*(ERImol(iorb2,iorb1,iorb1,iorb2)+ERImol(iorb1,iorb2,iorb2,iorb1))&
                   & -(ERImol(iorb2,iorb1,iorb2,iorb1)+ERImol(iorb1,iorb2,iorb1,iorb2)))
        enddo
        ! Inter PNOF7 ->  - sqrt(np hp nq hq) Lpq [Omegag and Omegaf, below and above Fl]
        if(ipnof==7) then
         do iorb2=nfl+1,forb-1
          horb2=1.0d0-OCC(iorb2)
          calc_epnof=calc_epnof-sqrt(OCC(iorb1)*OCC(iorb2)*horb1*horb2)&
                    & *(ERImol(iorb1,iorb2,iorb1,iorb2)+ERImol(iorb2,iorb1,iorb2,iorb1))
         enddo
         do iorb2=lorb+1,nbf
          horb2=1.0d0-OCC(iorb2)
          calc_epnof=calc_epnof-sqrt(OCC(iorb1)*OCC(iorb2)*horb1*horb2)&
                    & *(ERImol(iorb1,iorb2,iorb1,iorb2)+ERImol(iorb2,iorb1,iorb2,iorb1))
         enddo
        endif
       else
       ! SD energy ->   Jpp 
       calc_epnof=calc_epnof+ERImol(iorb1,iorb1,iorb1,iorb1)
       ! SD energy -> 1 nq (2 Jpq - Kpq)
       do iorb2=iorb1+1,nbf
        calc_epnof=calc_epnof+OCC(iorb2)&
                  & *(2.0d0*(ERImol(iorb2,iorb1,iorb1,iorb2)+ERImol(iorb1,iorb2,iorb2,iorb1))&
                  & -(ERImol(iorb2,iorb1,iorb2,iorb1)+ERImol(iorb1,iorb2,iorb1,iorb2)))
        enddo
       endif 
      enddo
      return
      end function
