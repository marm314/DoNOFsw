      SUBROUTINE EXTERN_OPT(COEF)
      USE PARCOM
      use m_noft_driver
      INTEGER::itermax=1000
      DOUBLE PRECISION::tolE=1.0d-8,tol_dif_Lambda=1.0d-4
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::COEF
      EXTERNAL::mo_ints
C-----------------------------------------------------------------------
      call run_noft(IPNOF,Ista,NBF,NBF5,NO1,NDOC,NCWO,NB,NA,ICGMETHOD,
     &  1,itermax,1,NTHRESHL,tolE,EN,COEF,mo_ints) 
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE EXTERN_OPT

      SUBROUTINE mo_ints(COEF,HCORE,ERImol)
      USE PARCOM
      USE PARCOM2
      DOUBLE PRECISION,DIMENSION(NBF,NBF),intent(inout)::HCORE
      DOUBLE PRECISION,DIMENSION(NBF,NBF,NBF,NBF),intent(inout)::ERImol
      DOUBLE PRECISION,DIMENSION(NBF,NBF),intent(in)::COEF
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      HCORE=0.0d0; ERImol=0.0d0;
      CALL HCOREc(COEF,HCORE)
      CALL ERIC1c(ERImol,IJKL,XIJKL,COEF,NBF)
      CALL ERIC23c(ERImol,COEF,NBF)
      CALL ERIC4c(ERImol,COEF,NBF)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE mo_ints

      SUBROUTINE HCOREc(COEF,HCORE)
      USE PARCOM
      USE PARCOM2
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::COEF,HCORE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::TEMPM
C-----------------------------------------------------------------------
      ALLOCATE(TEMPM(NBF,NBF))
      TEMPM=matmul(AHCORE2,COEF)
      HCORE=matmul(transpose(COEF),TEMPM)
      DEALLOCATE(TEMPM)
C-----------------------------------------------------------------------
      END SUBROUTINE HCOREc
