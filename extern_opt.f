      SUBROUTINE EXTERN_OPT(COEF,AHCORE,IERI,ERI)
      USE PARCOM
      USE PARCOM2
      use m_noft_driver
      INTEGER::i,j
      INTEGER,DIMENSION(NIJKL)::IERI
      DOUBLE PRECISION,DIMENSION(NIJKL)::ERI
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AHCORE
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::COEF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:)::ERImol
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::HCORE,TEMPM
C-----------------------------------------------------------------------
      ALLOCATE(HCORE(NBF,NBF),ERImol(NBF,NBF,NBF,NBF),TEMPM(NBF,NBF))
      AHCORE2=AHCORE
      TEMPM=0.0d0
      do j=1,NBF
       do i=1,NBF
        TEMPM(i,j)=COEF(i,j)
       enddo
      enddo
      CALL HCOREc(COEF,HCORE)
      CALL ERIC1c(ERImol,IERI,ERI,TEMPM,NBF)
      CALL ERIC23c(ERImol,TEMPM,NBF)
      CALL ERIC4c(ERImol,TEMPM,NBF)
      DEALLOCATE(TEMPM)
      if(ICGMETHOD==1) then
       call run_noft(IPNOF,Ista,NBF5,NO1,NDOC,NCWO,NB,NA,0,EN,
     & HCORE,ERImol) 
      else
       call run_noft(IPNOF,Ista,NBF5,NO1,NDOC,NCWO,NB,NA,1,EN,
     & HCORE,ERImol) 
       endif
      DEALLOCATE(HCORE,ERImol)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE EXTERN_OPT

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
