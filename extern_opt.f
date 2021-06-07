      SUBROUTINE EXTERN_OPT(COEF)
      USE M_VARS
      USE PARCOM
      USE PARCOM2
      use m_noft_driver
      LOGICAL::LRESTART
      INTEGER::itermax=1000,ie=0,nprt=0
      real(dp)::tolE=1.0d-8,tol_dif_Lambda=1.0d-4
      real(dp),DIMENSION(NBF,NBF)::COEF
      real(dp)::Enof
      real(dp),allocatable,dimension(:,:)::occup
      character(len=100)::ofile_name
      EXTERNAL::mo_ints
C-----------------------------------------------------------------------
      write(*,'(a)') 'Entering external RUN_NOFT optimization'
      allocate(occup(NBF,1))
      occup=0.0d0
      ofile_name='res.noft'
      call run_noft(INOF,Ista,NBF,NBF5,NO1,NDOC,NCWO,NB,NA,ie,ICGMETHOD,
     &  1,itermax,nprt,nprt,NTHRESHL,NDIIS,Enof,tolE,EN,COEF,OVERLAP2,
     &  occup(:,1),mo_ints,ofile_name)
!     &  occup(:,1),mo_ints,ofile_name,restart=LRESTART,ireadGAMMAS=1,
!     &  ireadOCC=1,ireadCOEF=1,ireadFdiag=1)
      deallocate(occup)
      write(*,'(a,f12.6)') 'Optimized energy',Enof
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE EXTERN_OPT

      SUBROUTINE mo_ints(NBFin,NBFin2,COEF,HCORE,ERImol)
      USE M_VARS
      USE PARCOM
      USE PARCOM2
      INTEGER,intent(in)::NBFin,NBFin2
      real(dp),DIMENSION(NBFin,NBFin),intent(inout)::HCORE
      real(dp),DIMENSION(NBFin,NBFin,NBFin,NBFin),intent(inout)::ERImol
      real(dp),DIMENSION(NBFin,NBFin),intent(in)::COEF
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      HCORE=0.0d0; ERImol=0.0d0;
      CALL HCOREc(COEF,HCORE)
      CALL ERIC1c(ERImol,IJKL,XIJKL,COEF,NBFin)
      CALL ERIC23c(ERImol,COEF,NBFin)
      CALL ERIC4c(ERImol,COEF,NBFin)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE mo_ints

      SUBROUTINE HCOREc(COEF,HCORE)
      USE M_VARS
      USE PARCOM
      USE PARCOM2
      real(dp),DIMENSION(NBF,NBF)::COEF,HCORE
      real(dp),ALLOCATABLE,DIMENSION(:,:)::TEMPM
C-----------------------------------------------------------------------
      ALLOCATE(TEMPM(NBF,NBF))
      TEMPM=matmul(AHCORE2,COEF)
      HCORE=matmul(transpose(COEF),TEMPM)
      DEALLOCATE(TEMPM)
C-----------------------------------------------------------------------
      END SUBROUTINE HCOREc

      SUBROUTINE elag_extern(COEF,RO,CJ12,CK12,ELAG)
      USE M_VARS
      USE PARCOM
      USE PARCOM2
      real(dp),DIMENSION(3)::DIPN
      real(dp),DIMENSION(NBF5)::RO
      real(dp),DIMENSION(NBF,NBF)::COEF,ELAG
      real(dp),DIMENSION(NBF5,NBF5)::CJ12,CK12
      real(dp),ALLOCATABLE,DIMENSION(:,:)::G,nCK12
      real(dp),ALLOCATABLE,DIMENSION(:,:,:)::QD
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE (G(NBF,NBF5),QD(NBF,NBF,NBF),nCK12(NBF5,NBF5))
      nCK12=-CK12
      CALL ENERGY1r(AHCORE2,IJKL,XIJKL,QD,COEF,RO,CJ12,nCK12,ELAG,
     &              DIPN,ADIPx,ADIPy,ADIPz,G)
      DEALLOCATE(G,QD)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE elag_extern
