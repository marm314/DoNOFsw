      SUBROUTINE EXTERN_OPT(COEF)
      USE PARCOM
      use m_noft_driver
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::COEF
      EXTERNAL::mo_ints1
C-----------------------------------------------------------------------
      if(ICGMETHOD==1) then
       write(*,*) 'Calling external module with CG'
       call run_noft(IPNOF,Ista,NBF,NBF5,NO1,NDOC,NCWO,NB,NA,0,EN,COEF,
     &    mo_ints1) 
      else
       write(*,*) 'Calling external module with LBFGS'
       call run_noft(IPNOF,Ista,NBF,NBF5,NO1,NDOC,NCWO,NB,NA,1,EN,COEF,
     &    mo_ints1) 
      endif
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE EXTERN_OPT

      SUBROUTINE mo_ints1(COEF,HCOREpp,ERI_J,ERI_K)
      USE PARCOM
      USE PARCOM2
      INTEGER::i,j,k
      DOUBLE PRECISION,DIMENSION(NBF5),intent(inout)::HCOREpp
      DOUBLE PRECISION,DIMENSION(NBFT5),intent(inout)::ERI_J,ERI_K
      DOUBLE PRECISION,DIMENSION(NBF,NBF),intent(in)::COEF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:)::ERImol
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::HCORE
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(HCORE(NBF,NBF),ERImol(NBF,NBF,NBF,NBF))
      HCORE=0.0d0; ERImol=0.0d0;
      CALL HCOREc(COEF,HCORE)
      CALL ERIC1c(ERImol,IJKL,XIJKL,COEF,NBF)
      CALL ERIC23c(ERImol,COEF,NBF)
      CALL ERIC4c(ERImol,COEF,NBF)
      k=1
      do i=1,NBF5
       HCOREpp(i)=HCORE(i,i)
       do j=1,i
        ERI_J(k)=ERImol(i,j,j,i) !J
        ERI_K(k)=ERImol(i,j,i,j) !K
        k=k+1
       enddo
      enddo
      DEALLOCATE(HCORE,ERImol)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE mo_ints1

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
