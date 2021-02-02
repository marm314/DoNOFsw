!======================================================================!
!                                                                      !
!              N O F - M B P T  S U B R O U T I N E S                  !
!                                                                      !
!                                                                      !
!======================================================================!
! ==================================================================== !
!                           Date: 10/01/2021                           !
!                                                                      !
! mbpt.f contains the subroutines related to the RPA+GW+SOSEX          ! 
! calculations. RPA is used for OMEGAs and Sigma_c for Ec              !
!                                                                      !
!  Fpq: Form Fpq factor from SOSEX see X. Ren. PRB 88, 035120 (2013)   !
!  Wpqrs: Screened Coulomb, see MOLGW paper                            !
!  Integrated_omega: See RPA paper, X. PRB 88, 035120 (2013)           !
!                                                                      !
! ==================================================================== !
C MBPT
      SUBROUTINE MBPTCALC(ELAG,COEF,RO,CJ12,CK12,AHCORE,ADIPx,ADIPy,
     & ADIPz,IERI,ERI)
      USE PARCOM
      PARAMETER (tol10=1.0D-10)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (FOURTH=0.25D0)
      PARAMETER (ONE=1.0D0)
      PARAMETER (TWO=2.0D0)
      PARAMETER (FOUR=4.0D0)
      LOGICAL::diagFOCK
      INTEGER::i,j,k,l,a,b,info,last_coup
      INTEGER::order
      INTEGER,DIMENSION(NIJKL)::IERI
      DOUBLE PRECISION::AUtoEV,EPNOF,ESDc,EcRPA,EcMP2,EcGoWo,EcGMSOS
      DOUBLE PRECISION::EcGoWoSOS 
      DOUBLE PRECISION::iEcRPA,iEcSOSEX,iEcRPASOS,EcCCSD,mu
      DOUBLE PRECISION,DIMENSION(NIJKL)::ERI
      DOUBLE PRECISION,DIMENSION(NBF5)::RO
      DOUBLE PRECISION,DIMENSION(NBF5,NBF5)::CJ12,CK12
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::AHCORE,ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NBF,NBF)::ELAG,COEF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::OCC,EIG,BIGOMEGA,TEMPV
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::CINTER,CINTRA 
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::TEMPM,TEMPM2,FOCKm
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::ApB,AmB
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::XmY,XpY 
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:)::ERImol,ERImol2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:)::wmn,wmn2
C-----------------------------------------------------------------------
C     NCO:  Number of HF occupied MOs (OCC=1 in SD)
C     NVIR: Number of HF virtual  MOs (OCC=0 in SD)
C
C     NO1:  Number of inactive doubly occupied orbitals (OCC=1)
C     NDOC: Number of strongly occupied MOs
C     NCWO: Number of coupled weakly occupied MOs per strongly occupied
C     NCWO*NDOC: Active orbitals in the virtual subspace
C     NO0: Empty orbitals  (OCC=0)
C     NAC:  Dimension of the active natural orbital subspace
C
C     NB/NCO: Number of orbitals in the doubly-occupied space
C     NA: NB/NCO + Number of orbitals in the singly-occupied space (NSOC)
C
C           NCO + NSOC     |       NVIR          = NBF
C-----------------------------------------------------------------------
      IF(ICOEF==0)CALL ELAGCOEF0(ELAG,COEF,RO,CJ12,CK12,AHCORE,
     &                           ADIPx,ADIPy,ADIPz,IERI,ERI)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Initial allocation and def. of number of pairs (occ,virtual) i <-> a  
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      AUtoEV=27.211399
      NVIR=NBF-NA
      Nab=NVIR*NA !NVIR*NCO
      last_coup=NA+NCWO*(NCO-NO1PT2)
      WRITE(6,1)NBF,NO1PT2,NA,NVIR,last_coup,Nab
      ALLOCATE(OCC(NBF),TEMPM(NBF,NBF))
      OCC=ZERO
      DO i=1,NBF5
       OCC(i) = RO(i)
      ENDDO
      DO j=1,NBF
       do i=1,NBF
        TEMPM(i,j) = COEF(i,j)
       enddo
      ENDDO
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Prepare coef. factors CINTRA and CINTER.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(TUNEMBPT) then ! as used in PRL and PRA papers for NCWO=-1
       ALLOCATE(CINTER(NBF),CINTRA(NBF))
       call tune_coefs(CINTER,CINTRA,OCC,NBF,NA)
      endif
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  One orbital energies (EIG) and FOCK matrix (FOCKm) with COEF/=CHF
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(EIG(NBF),FOCKm(NBF,NBF))
      CALL FOCKMOL(NBF,COEF,TEMPM,ELAG,EIG,FOCKm,AHCORE,IERI,ERI,ESD)
      diagFOCK=.true.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Tune Fock(i,a) elements and calc. 
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(TUNEMBPT) then
       ALLOCATE(TEMPM2(NBF,NBF))
       call tunefock(CINTRA,CINTER,FOCKm,NBF,NCO,NVIR,NCWO,NO1PT2,EIG,
     &  TEMPM2,diagFOCK,NA)
      endif
      DEALLOCATE(FOCKm)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Print orb. energies and occ numbers.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      write(*,*) '          '
      write(*,*) ' List of qp-orbital energies (a.u.) and occ numbers 
     & used'
      write(*,*) '          '
      mu=(EIG(NA+1)+EIG(NA))/TWO
      DO i=1,NBF
       EIG(i)=EIG(i)-mu
       write(6,2) EIG(i),OCC(i)
      ENDDO
      write(*,*) '          '
      write(*,'(a,f10.5)') ' Chemical potential used for qp-orbital 
     & energies',mu
      write(*,*) '          '
      DEALLOCATE(OCC) ! We do not need the OCCs anymore
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Allocate 2e_integrals array(s) and form ERI in MO basis (see mp2.f)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ALLOCATE(ERImol(NBF,NBF,NBF,NBF))
      IF(TUNEMBPT) ALLOCATE(ERImol2(NBF,NBF,NBF,NBF))
      CALL ERIC1c(ERImol,IERI,ERI,TEMPM,NBF)
      CALL ERIC23c(ERImol,TEMPM,NBF)
      CALL ERIC4c(ERImol,TEMPM,NBF)
      IF(TUNEMBPT) ERImol2=ERImol
      DEALLOCATE(TEMPM)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Print 2e- integrals in MO basis
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     call print2eint(NBF,ERImol)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Modify ERImol with CINTER and CINTRA before computing W, wmn, etc.
C  Also transform ERIs from NO basis to "DIAG(FOCK)" basis (TEMPM2)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(TUNEMBPT) then   
       call tuneerimol(CINTER,CINTRA,NBF,NCO,NVIR,NCWO,NO1PT2,ERImol,
     &  TDHF,NA)
       DEALLOCATE(CINTER,CINTRA)
       if(diagFOCK.eqv..false.) then
        if(mbptmem) then
         call transform2mem(NBF,TEMPM2,ERImol)
         call transform2mem(NBF,TEMPM2,ERImol2)
        else
         call transform2(NBF,TEMPM2,ERImol,ERImol2)
        endif
       endif
       DEALLOCATE(TEMPM2)
      endif
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Build A+B and A-B matrices for CASIDA EQN and 1st contribution to EcRPA.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      EcRPA=ZERO
      ALLOCATE(AmB(Nab,Nab),ApB(Nab,Nab))
      AmB=ZERO
      ApB=ZERO
      k=1      
      do i=1,NA
        do a=NA+1,nbf
          l=1
          do j=1,NA
           do b=NA+1,nbf
             ApB(k,l)=FOUR*ERImol(i,j,b,a) ! This is i-a and j-b, the 4 is from the spin sum 
             AmB(k,l)=ZERO 
            if(i.eq.j.and.a.eq.b) then      
             AmB(k,l)=EIG(a)-EIG(i)
             ApB(k,l)=ApB(k,l)+AmB(k,l)
            endif 
            if(TDHF) then
             ApB(k,l)=ApB(k,l)-ERImol(i,a,b,j)-ERImol(i,a,j,b) ! i-j,a-b and i-b,a-j
             AmB(k,l)=AmB(k,l)-ERImol(i,a,b,j)+ERImol(i,a,j,b) ! i-j,a-b and i-b,a-j
            endif
            if(k.eq.l) then ! 1st contrib. to RPA -> 0.25*diag_terms  
             EcRPA=EcRPA-FOURTH*(ApB(k,k)+AmB(k,k)) 
            endif  
            l=l+1
           enddo
          enddo
         k=k+1
        enddo
      enddo
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Solve Casida Equation by Cholesky Decomposition (calling LAPACK):
C  sqrt(BIGOMEGA)=Excitation energies
C  the eigenvectors (F) of C*F=BIGOMEGA^2 *F are stored in columns in TEMPM
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      info=0
      CALL DPOTRF('L',Nab,ApB,Nab,info)
      if(info.ne.0) then
       write(*,*) ' ' 
       write(*,*) ' WARNING! DPOTRF FAILED TO CONVERGE ',info 
       write(*,*) ' '  
      else
       CALL DSYGST(3,'L',Nab,AmB,Nab,ApB,Nab,info)
       if(info.ne.0) then 
        write(*,*) ' ' 
        write(*,*) ' WARNING! DSYGST FAILED TO CONVERGE ',info
        write(*,*) ' '
       endif
      endif
      ALLOCATE(TEMPM(Nab,Nab),BIGOMEGA(Nab),TEMPV(Nab))
      ALLOCATE(XpY(Nab,Nab),XmY(Nab,Nab))
      do i=1,Nab
       do j=i+1,Nab
         AmB(i,j)=AmB(j,i)
       enddo
      enddo
      BIGOMEGA=ZERO
      TEMPV=ZERO 
      CALL DIAG(Nab,AmB,TEMPM,BIGOMEGA,TEMPV) 
      BIGOMEGA(:) = SQRT(BIGOMEGA(:))
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Compute eigenvectors X+Y and X-Y as in MOLGW (calling LAPACK)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      XpY=TEMPM
      XmY=TEMPM
      if(info.eq.0) then
       call DTRMM('L','L','N','N',Nab,Nab,ONE,ApB,Nab,XmY,Nab)
       call DTRSM('L','L','T','N',Nab,Nab,ONE,ApB,Nab,XpY,Nab)
      endif
      forall(i=1:Nab)
       XpY(:,i)=XpY(:,i)*SQRT(BIGOMEGA(i))
       XmY(:,i)=XmY(:,i)/SQRT(BIGOMEGA(i))
      endforall
      DEALLOCATE(AmB,ApB,TEMPM,TEMPV,XmY) ! We do not need X-Y anymore 
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Compute omegas, oscillator strenghts and static polarizability (omega->0)       
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call td_polarizability(NBF,NCO,Nab,NA,COEF,XpY,BIGOMEGA,
     & ADIPx,ADIPy,ADIPz,EcRPA,TDHF,AUtoEV)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GoWo Ec energy using Galitskii-Migdal EQN
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Compute w^s mn = sum_ia <im|an>* (X^s _ia + Y^s_ia) -> see MolGW paper (i is occ, a is virt) 
      ALLOCATE(wmn(NBF,NBF,Nab))
      IF(TUNEMBPT) ALLOCATE(wmn2(NBF,NBF,Nab))
      write(*,*) ' '
      if(TUNEMBPT) then   ! Do GW for NOFc-GW
       write(*,*) 'Computing GW correction for NOFc-GW'
      else                        ! Do standard GW NOT FOR NOF
       write(*,*) 'Computing GW standard correction'
      endif
      wmn=ZERO
      IF(TUNEMBPT) wmn2=ZERO
      call build_wmn(NBF,Nab,NA,wmn,wmn2,TUNEMBPT,ERImol,ERImol2,XpY)
      EcGoWo=ZERO
      EcGMSOS=ZERO
      if(TUNEMBPT) then
      call gw_gm_eq(NA,NCO,NBF,Nab,wmn,wmn2,EIG,EcGoWo,EcGMSOS,
     &  TUNEMBPT,XpY,BIGOMEGA,ERImol2)
      else
      call gw_gm_eq(NA,NCO,NBF,Nab,wmn,wmn2,EIG,EcGoWo,EcGMSOS,
     &  TUNEMBPT,XpY,BIGOMEGA,ERImol)
      endif
      EcGoWo=TWO*EcGoWo
      EcGoWoSOS=EcGoWo+EcGMSOS
      DEALLOCATE(XpY)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  MP2 Ec energy 
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(TUNEMBPT) then
       write(*,*) 'Computing MP2 correction for NOFc-MP2'
      else
       write(*,*) 'Computing MP2 standard correction'
      endif
      EcMP2=ZERO
      call mp2_eq(NA,NCO,NBF,EIG,ERImol,ERImol2,TUNEMBPT,EcMP2)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Integrated Ec RPA, AC-SOSEX, and RPA+AC-SOSEX, etc. 
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(TUNEMBPT) then
       write(*,*) 'Computing integrated RPA+AC-SOSEX correction 
     &for NOFc-RPA and NOFc-RPA+AC-SOSEX'
      else
       write(*,*)'Computing integrated RPA+AC-SOSEX standard correction'
      endif
      iEcRPA=0.0d0
      iEcSOSEX=0.0d0
      call rpa_sosex_eq(NA,NCO,NBF,Nab,ERImol,ERImol2,EIG,BIGOMEGA,wmn,
     &  iEcRPA,iEcSOSEX,TUNEMBPT)
      iEcRPASOS=iEcRPA+iEcSOSEX
      DEALLOCATE(wmn,BIGOMEGA)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Ec CCSD 
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(TUNEMBPT) then
       write(*,*) 'Computing CCSD correction for NOFc-CCSD'
      else
       write(*,*)'Computing CCSD standard correction'
      endif
      EcCCSD=0.0d0
      call ccsd_eq(NA,NCO,NBF,ERImol,ERImol2,EIG,EcCCSD,TUNEMBPT,
     &  CCSD_READ)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Write  Ec energies and final results
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ESD=ESD+EN+ECndHF
      if(TUNEMBPT) then
       ESDc=ESD+ECndl
      else
       ESDc=ESD
      endif
      EPNOF=EELEC+EN
      IF(TUNEMBPT) THEN ! SD + DYN + ND (can be used with CANONICAL orbs with ICOEF=0)
      write(6,3)ESD,ESDc,EPNOF,ECndl,EcRPA,iEcRPA,
     & iEcSOSEX,iEcRPASOS,EcGoWo,EcGMSOS,EcGoWoSOS,EcMP2,EcCCSD,
     & ESD+EcRPA,ESD+iEcRPA,
     & ESD+iEcRPASOS,ESD+EcGoWo,ESD+EcGoWoSOS,
     & ESD+EcMP2,ESD+EcCCSD,ESDc+EcRPA,ESDc+iEcRPA,ESDc+iEcRPASOS,
     & ESDc+EcGoWo,ESDc+EcGoWoSOS,ESDc+EcMP2,ESDc+EcCCSD
      ENDIF
      IF(TUNEMBPT.eqv..false.) then ! HF(NOF orbs) + DYN (NOF orbs = CANONICAL orbs for ICOEF=0)
      write(6,4)ESD,EPNOF,EcRPA,iEcRPA,iEcSOSEX,iEcRPASOS,
     & EcGoWo,EcGMSOS,EcGoWoSOS,EcMP2,EcCCSD,ESD+EcRPA,ESD+iEcRPA,
     & ESD+iEcRPASOS,ESD+EcGoWo,ESD+EcGoWoSOS,ESD+EcMP2,ESD+EcCCSD
      ENDIF
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  Deallocate and clean 
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(TUNEMBPT) DEALLOCATE(ERImol2,wmn2)
      DEALLOCATE(EIG,ERImol)
      RETURN
C-----------------------------------------------------------------------
    1 FORMAT(2X,//,3X,'MBPT ENERGIES',/
     *       2X,15('='),//,
     *       1X,'Number of orbitals        (NBASIS) =',I5,/,
     *       1X,'Number of frozen pairs       (NFR) =',I5,/,
     *       1X,'Number of occupied orbs.     (NOC) =',I5,/,
     *       1X,'Number of virtual orbs.     (NVIR) =',I5,/, 
     *       1X,'Last coupled orbital        (NLAS) =',I5,/, 
     *       1X,'Size of A+B and A-B (NAB=NOCxNVIR) =',I5) 
    2 FORMAT(3X,F15.10,' ',F15.10,' ')
    3 FORMAT(/,1X,'E(SD)                  ',5X,F20.10,' a.u.',/,
     * 1X,'E(SD+ND)               ',5X,F20.10,' a.u.',/,
     * 1X,'E(PNOFi)               ',5X,F20.10,' a.u.',/,
     * ' ',/,
     * 1X,'Ec(ND)                 ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(RPA-FURCHE)         ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(RPA)                ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(AC-SOSEX)           ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(RPA+AC-SOSEX)       ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(GW@GM)              ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(SOSEX@GM)           ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(GW@GM+SOSEX@GM)     ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(MP2)                ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(CCSD)               ',5X,F20.10,' a.u.',/,
     * ' ',/,
     * 1X,'E(RPA-FURCHE)          ',5X,F20.10,' a.u.',/,
     * 1X,'E(RPA)                 ',5X,F20.10,' a.u.',/,
     * 1X,'E(RPA+AC-SOSEX)        ',5X,F20.10,' a.u.',/,
     * 1X,'E(GW@GM)               ',5X,F20.10,' a.u.',/,
     * 1X,'E(GW@GM+SOSEX@GM)      ',5X,F20.10,' a.u.',/,
     * 1X,'E(MP2)                 ',5X,F20.10,' a.u.',/,
     * 1X,'E(CCSD)                ',5X,F20.10,' a.u.',/,
     * ' ',/,
     * 1X,'E(NOFc-RPA-FURCHE)     ',5X,F20.10,' a.u.',/,
     * 1X,'E(NOFc-RPA)            ',5X,F20.10,' a.u.',/,
     * 1X,'E(NOFc-RPA+AC-SOSEX)   ',5X,F20.10,' a.u.',/,
     * 1X,'E(NOFc-GW@GM)          ',5X,F20.10,' a.u.',/,
     * 1X,'E(NOFc-GW@GM+SOSEX@GM) ',5X,F20.10,' a.u.',/,
     * 1X,'E(NOFc-MP2)            ',5X,F20.10,' a.u.',/,
     * 1X,'E(NOFc-CCSD)           ',5X,F20.10,' a.u.')
    4 FORMAT(/,1X,'E(SD)                 ',5X,F20.10,' a.u.',/,
     * 1X,'E(PNOFi)              ',5X,F20.10,' a.u.',/,
     * ' ',/,
     * 1X,'Ec(RPA-FURCHE)        ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(RPA)               ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(AC-SOSEX)          ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(RPA+AC-SOSEX)      ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(GW@GM)             ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(SOSEX@GM)          ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(GW@GM+SOSEX@GM)    ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(MP2)               ',5X,F20.10,' a.u.',/,
     * 1X,'Ec(CCSD)              ',5X,F20.10,' a.u.',/,
     * ' ',/,
     * 1X,'E(RPA-FURCHE)         ',5X,F20.10,' a.u.',/,
     * 1X,'E(RPA)                ',5X,F20.10,' a.u.',/,
     * 1X,'E(RPA+AC-SOSEX)       ',5X,F20.10,' a.u.',/,
     * 1X,'E(GW@GM)              ',5X,F20.10,' a.u.',/,
     * 1X,'E(GW@GM+SOSEX@GM)     ',5X,F20.10,' a.u.',/,
     * 1X,'E(MP2)                ',5X,F20.10,' a.u.',/,
     * 1X,'E(CCSD)               ',5X,F20.10,' a.u.')
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE MBPTCALC

      function integrated_omega(i,j,a,b,order,nbf,nab,
     & wmn,weights,freqs,cfreqs,EIG,BIGOMEGA,ERImolijab)
      implicit none
      integer::ii
      integer,intent(in)::i,j,a,b,nbf,nab,order
      double precision::integrated_omega,Fpq
      double precision,dimension(nbf),intent(in)::EIG
      double precision,dimension(nab)::BIGOMEGA
      double precision,dimension(nbf,nbf,nab),intent(in)::wmn
      double precision,dimension(order),intent(in)::freqs,weights
      double precision,intent(in)::ERImolijab
      complex,dimension(order),intent(in)::cfreqs
      double precision::tol4
      complex::integral,wpqrsc
      tol4=1.0D-4
      integral=0.0d0
      integrated_omega=0.0d0
      do ii=1,order
       call Wpqrs(i,j,a,b,nab,nbf,wmn,cfreqs(ii),BIGOMEGA,ERImolijab,
     & wpqrsc)   
       integral=integral+weights(ii)*
     & Fpq(i,a,nbf,EIG,freqs(ii))*Fpq(j,b,nbf,EIG,freqs(ii))*
     & wpqrsc
      enddo
      integrated_omega=REAL(integral)
      if(AIMAG(integral).gt.tol4) then
       write(6,12) i,j,a,b,AIMAG(integral)
      endif
      return
   12 FORMAT('WARNING! LARGE IMAGINARY FREQ. INTEG. IN',1X,I5,1X,I5,
     & 1X,I5,1X,I5,5X,F15.10)
      end function
 
      function Fpq(p,q,nbf,EIG,freq)
      implicit none
      integer,intent(in)::p,q,nbf
      double precision,intent(in)::freq
      double precision,dimension(nbf),intent(in)::EIG
      double precision::Fpq
       Fpq=2.0d0*(EIG(p)-EIG(q))/((EIG(p)-EIG(q))**2.0d0+freq**2.0d0)
      return
      end function
     
      subroutine Wpqrs(p,q,r,s,nab,nbf,wmn,freqc,BIGOMEGA,ERImolpqrs,
     &  wpqrsc)
      implicit none
      integer::ii
      integer,intent(in)::p,q,r,s,nab,nbf
      double precision,intent(in)::ERImolpqrs
      double precision,dimension(nab),intent(in)::BIGOMEGA
      double precision,dimension(nbf,nbf,nab),intent(in)::wmn
      complex,intent(in)::freqc
      complex::wpqrsc
      wpqrsc=ERImolpqrs 
      do ii=1,nab                      
       wpqrsc=wpqrsc+wmn(p,r,ii)*wmn(q,s,ii)*(1.0d0/(freqc-BIGOMEGA(ii))
     & -1.0d0/(freqc+BIGOMEGA(ii)))
      enddo
      end subroutine Wpqrs
      
      subroutine tune_coefs(CINTER,CINTRA,OCC,NBF,NA)
      implicit none
      integer,intent(in)::NBF,NA
      double precision,dimension(nbf),intent(in)::OCC
      double precision,dimension(nbf),intent(inout)::CINTER,CINTRA
      integer::ii
      double precision::hi
      CINTER=1.0d0
      CINTRA=1.0d0
      do ii=1,NBF
       hi=1.0d0-OCC(ii)
       ! CINTRA 
       CINTRA(ii)=1.0d0-DABS(1.0d0-2.0d0*OCC(ii))
       CINTRA(ii)=1.0d0-CINTRA(ii)**2.0d0
       ! CINTER
       if(ii.le.NA) then ! eq. (34) PRA multiplets
         CINTER(ii)=1.0d0
       else
        CINTER(ii)=DABS(1.0d0-2.0d0*OCC(ii))
        CINTER(ii)=CINTER(ii)**2.0d0
       endif
      enddo
      write(*,*) ' '
      end subroutine tune_coefs

      subroutine tunefock(CINTRA,CINTER,FOCKm,NBF,NCO,NVIR,NCWO,NFR,EIG,
     & EIGENVE,diagF,NA)
      implicit none
      logical,intent(inout)::diagF
      integer,intent(in)::NBF,NCO,NVIR,NCWO,NFR,NA
      double precision,dimension(NBF),intent(inout)::CINTER,CINTRA,EIG
      double precision,dimension(NBF,NBF),intent(inout)::FOCKm,EIGENVE
      integer,allocatable,dimension(:)::coup
      double precision,allocatable,dimension(:)::TEMPV
      double precision::tol6
      integer::i,k,in,kn,lmin_i,lmax_i,last_coup
      allocate(coup(NVIR))
      coup=0
      tol6=1.0D-6
      do i=1,NBF-1
       do k=i+1,NBF
        if(abs(FOCKm(i,k)).gt.tol6) then
         diagF=.false.
         goto 687
        endif
       enddo
      enddo
  687 continue
      if(diagF.eqv..false.) then
       last_coup=NA+ncwo*(nco-nfr)
       allocate(TEMPV(NBF))
       do i=1,nco
        lmin_i = NA+ncwo*(nco-i)+1
        lmax_i = NA+ncwo*(nco-i)+ncwo
        do k=1,nvir
         kn=k+NA
         FOCKm(i,kn)=0.0d0 
         FOCKm(kn,i)=0.0d0 
         if((lmin_i<=kn.and.kn<=lmax_i).and.(lmax_i<=last_coup)) then
          coup(k)=i !coup(kn)=i
         endif
        enddo
       enddo
       do i=nco+1,NA
        do k=1,nvir
         kn=k+NA
         FOCKm(i,kn)=0.0d0
         FOCKm(kn,i)=0.0d0
        enddo
       enddo
       do i=1,NA
        do k=i+1,NA
         FOCKm(i,k)=FOCKm(i,k)*CINTER(i)*CINTER(k)
         FOCKm(k,i)=FOCKm(i,k)
        enddo
       enddo
       do i=1,nvir
        in=i+NA
        do k=i+1,nvir
         kn=k+NA
         if(coup(i)==coup(k)) then
          FOCKm(in,kn)=FOCKm(in,kn)*CINTRA(in)*CINTRA(kn)
         else
          FOCKm(in,kn)=FOCKm(in,kn)*CINTER(in)*CINTER(kn)
         endif
         FOCKm(kn,in)=FOCKm(in,kn)
        enddo
       enddo 
       write(*,*) ' '
       CALL DIAG(NBF,FOCKm,EIGENVE,EIG,TEMPV) 
       deallocate(TEMPV)
      else
       EIGENVE=0.0d0
       do i=1,nbf
        EIGENVE(i,i)=1.0d0
       enddo
      endif
      deallocate(coup)
      end subroutine tunefock
     
      subroutine tuneerimol(CINTER,CINTRA,nbf,nco,nvir,ncwo,nfr,ERImol,
     &  TDHF,NA)
      implicit none
      logical,intent(in)::TDHF
      integer,intent(in)::NBF,NCO,NVIR,NCWO,NFR,NA
      double precision,dimension(NBF),intent(in)::CINTER,CINTRA
      double precision,dimension(nbf,nbf,nbf,nbf),intent(inout)::ERImol
      integer::i,j,k,l,kn,ln,lmin_i,lmax_i,last_coup
      double precision::Ciikl,Cijkl
      last_coup=NA+ncwo*(nco-nfr)
      do i=1,NA
       lmin_i = NA+ncwo*(nco-i)+1
       lmax_i = NA+ncwo*(nco-i)+ncwo
       do j=1,NA
        ! Subspace 'i' occ, still could be inter or intra
        if(j==i.and.j<=nco.and.i<=nco) then
         do k=1,nvir
          kn = k + NA
          do l=1,nvir
           ln = l + NA
           if(i>NCO .or. j>NCO)then
            Ciikl = CINTER(kn)*CINTER(ln)*CINTER(i)*CINTER(i)
           elseif((     (lmin_i<=kn.and.kn<=lmax_i)
     &        .and. (lmin_i<=ln.and.ln<=lmax_i)) 
     &        .and. (lmax_i<=last_coup)        )then
            Ciikl = CINTRA(kn)*CINTRA(ln)*CINTRA(i)*CINTRA(i)
           else
            Ciikl = CINTER(kn)*CINTER(ln)*CINTER(i)*CINTER(i)
           endif
           ERImol(i,i,ln,kn) = Ciikl*ERImol(i,i,ln,kn)
           ERImol(kn,i,ln,i) = ERImol(i,i,ln,kn)
           ERImol(kn,ln,i,i) = ERImol(i,i,ln,kn)
           ERImol(i,ln,i,kn) = ERImol(i,i,ln,kn)
           ERImol(i,i,kn,ln) = ERImol(i,i,ln,kn)
           ERImol(ln,i,kn,i) = ERImol(i,i,ln,kn)
           ERImol(ln,kn,i,i) = ERImol(i,i,ln,kn)
           ERImol(i,kn,i,ln) = ERImol(i,i,ln,kn)
           if(TDHF) then
            ERImol(i,kn,ln,i) = Ciikl*ERImol(i,kn,ln,i)
            ERImol(i,ln,kn,i) = ERImol(i,kn,ln,i)
            ERImol(ln,i,i,kn) = ERImol(i,kn,ln,i)
            ERImol(kn,i,i,ln) = ERImol(i,kn,ln,i)
           endif
          enddo
         enddo
        else          ! Purely interspace
         do k=1,nvir
          kn = k + NA
          do l=1,nvir
           ln = l + NA
           Cijkl = CINTER(kn)*CINTER(ln)*CINTER(i)*CINTER(j)
           ERImol(i,j,ln,kn) = Cijkl*ERImol(i,j,ln,kn)
           ERImol(i,ln,j,kn) = ERImol(i,j,ln,kn)
           ERImol(kn,ln,j,i) = ERImol(i,j,ln,kn)
           ERImol(kn,j,ln,i) = ERImol(i,j,ln,kn)
           ERImol(j,i,kn,ln) = ERImol(i,j,ln,kn)
           ERImol(j,kn,i,ln) = ERImol(i,j,ln,kn)
           ERImol(ln,kn,i,j) = ERImol(i,j,ln,kn)
           ERImol(ln,i,kn,j) = ERImol(i,j,ln,kn)
           if(TDHF) then
            ERImol(i,kn,ln,j) = Cijkl*ERImol(i,kn,ln,j)
            ERImol(i,ln,kn,j) = ERImol(i,kn,ln,j)
            ERImol(j,ln,kn,i) = ERImol(i,kn,ln,j)
            ERImol(j,kn,ln,i) = ERImol(i,kn,ln,j)
            ERImol(kn,j,i,ln) = ERImol(i,kn,ln,j)
            ERImol(kn,i,j,ln) = ERImol(i,kn,ln,j)
            ERImol(ln,i,j,kn) = ERImol(i,kn,ln,j)
            ERImol(ln,j,i,kn) = ERImol(i,kn,ln,j)
           endif
          enddo
         enddo
        endif
       enddo
      enddo
      end subroutine tuneerimol      

      subroutine td_polarizability(NBF,NCO,Nab,NA,COEF,XpY,
     & BIGOMEGA,ADIPx,ADIPy,ADIPz,EcRPA,TDHF,AUtoEV)
      implicit none
      LOGICAL,intent(in)::TDHF
      INTEGER,intent(in)::Nbf,Nab,NCO,NA
      DOUBLE PRECISION,intent(inout)::EcRPA
      DOUBLE PRECISION,intent(in)::AUtoEV
      DOUBLE PRECISION,DIMENSION(NBF,NBF),intent(in)::ADIPx,ADIPy,ADIPz
      DOUBLE PRECISION,DIMENSION(NBF,NBF),intent(in)::COEF
      DOUBLE PRECISION,DIMENSION(Nab),intent(in)::BIGOMEGA
      DOUBLE PRECISION,DIMENSION(Nab,Nab),intent(in)::XpY
      INTEGER::i,j,k,l,a
      DOUBLE PRECISION::tol10,tol6
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) ::MDIPx,MDIPy,MDIPz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) ::COEFT,TEMPM
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) ::STATICPOL
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) ::OSCSTR,DIPSUM
      tol6=1.0D-6
      tol10=1.0D-10
C  Prepare dipole moment from AO to MO
      ALLOCATE(MDIPx(NBF,NBF),MDIPy(NBF,NBF),MDIPz(NBF,NBF))
      ALLOCATE(COEFT(NBF,NBF),TEMPM(NBF,NBF))
      do i=1,NBF
       do j=1,NBF
        COEFT(i,j)=COEF(j,i)  ! MO COEFs in columns
       enddo
      enddo
C      x
      TEMPM=matmul(ADIPx,COEF) 
      TEMPM=matmul(COEFT,TEMPM)
      MDIPx=TEMPM
C      y
      TEMPM=matmul(ADIPy,COEF) 
      TEMPM=matmul(COEFT,TEMPM)
      MDIPy=TEMPM
C      z
      TEMPM=matmul(ADIPz,COEF) 
      TEMPM=matmul(COEFT,TEMPM)
      MDIPz=TEMPM
      DEALLOCATE(COEFT,TEMPM)
C  Compute OSCILLATOR STRENGHTS and 2nd contribution to EcRPA
      ALLOCATE(OSCSTR(Nab),DIPSUM(3),TEMPM(3,Nab))
      OSCSTR=0.0d0
      do k=1,Nab
       EcRPA=EcRPA+0.5d0*ABS(BIGOMEGA(k))  ! 2nd contribution to EcRPA
       DIPSUM=0.d0
       l=1
       do i=1,NA
        do a=NA+1,NBF
          DIPSUM(1)=DIPSUM(1)+MDIPx(i,a)*XpY(l,k)*DSQRT(2.0d0)
          DIPSUM(2)=DIPSUM(2)+MDIPy(i,a)*XpY(l,k)*DSQRT(2.0d0)
          DIPSUM(3)=DIPSUM(3)+MDIPz(i,a)*XpY(l,k)*DSQRT(2.0d0)
          l=l+1
        enddo
       enddo
       TEMPM(1,k)=DIPSUM(1)
       TEMPM(2,k)=DIPSUM(2)
       TEMPM(3,k)=DIPSUM(3)
       OSCSTR(k)=2.0d0*(DIPSUM(1)**2.0d0+DIPSUM(2)**2.0d0+
     & DIPSUM(3)**2.0d0)*BIGOMEGA(k)/3.0d0
      enddo
C  Print Omegas and Oscillator strenghts
      write(*,*) ' '
      if(TDHF) then
       write(*,*) 'TD-HF CASIDA eq. solved'
      else 
       write(*,*) 'TD-H (RPA) CASIDA eq. solved' 
      endif
      write(*,*) 'N. excitation   a.u.         eV            nm      
     &osc. strenght'
      do i=1,Nab
       if(OSCSTR(i).gt.tol6) then
        write(6,13) i,BIGOMEGA(i),BIGOMEGA(i)*AUtoEV,1239.84193/
     &  (BIGOMEGA(i)*AUtoEV),OSCSTR(i)
        endif
      enddo
      DEALLOCATE(OSCSTR,MDIPx,MDIPy,MDIPz,DIPSUM)
C  Static polarizability alpha_xx' = <x|Xi(r,r',w=0)|x'> 
      ALLOCATE(STATICPOL(3,3))
      STATICPOL=0.0d0
      do i=1,Nab
       forall(j=1:3, k=1:3)
            STATICPOL(j,k)=STATICPOL(j,k) 
     &     + 2.0d0*TEMPM(j,i)*TEMPM(k,i)/(BIGOMEGA(i)+tol10)
       end forall
      enddo
      write(*,*) ' '
      write(6,14) STATICPOL(1,1),STATICPOL(1,2),STATICPOL(1,3), 
     &   STATICPOL(2,1),STATICPOL(2,2),STATICPOL(2,3),
     &   STATICPOL(3,1),STATICPOL(3,2),STATICPOL(3,3)
      write(6,15) ((STATICPOL(1,1)+STATICPOL(2,2)+STATICPOL(3,3))/3.0d0)
      DEALLOCATE(TEMPM,STATICPOL)

   13 FORMAT(3X,I5,3X,F12.5,1X,F12.5,1X,F12.4,1X,F12.6)
   14 FORMAT(1X,'Static polarizability:',//,
     *       F15.10,3X,F15.10,3X,F15.10,/,
     *       F15.10,3X,F15.10,3X,F15.10,/,
     *       F15.10,3X,F15.10,3X,F15.10)
   15 FORMAT(/,1X,'Trace of the static polarizability',F15.10)

      end subroutine td_polarizability

      subroutine build_wmn(NBF,Nab,NA,wmn,wmn2,TUNEMBPT,ERImol,ERImol2,
     &  XpY)
      implicit none
      logical,intent(in)::TUNEMBPT
      integer,intent(in)::NBF,Nab,NA
      double precision,dimension(NBF,NBF,Nab),intent(inout)::wmn,wmn2 
      double precision,dimension(Nab,Nab),intent(in)::XpY
      double precision,dimension(NBF,NBF,NBF,NBF),intent(in)::ERImol
      double precision,dimension(NBF,NBF,NBF,NBF),intent(in)::ERImol2
      integer::i,j,k,l,a,b
      do k=1,Nab
       do i=1,NBF
        do j=1,i
         a=1
         b=NA+1
         do l=1,Nab
          wmn(i,j,k)=wmn(i,j,k)+ERImol(a,i,j,b)*XpY(l,k)
          IF(TUNEMBPT) THEN
           wmn2(i,j,k)=wmn2(i,j,k)+ERImol2(a,i,j,b)*XpY(l,k)
          ENDIF
          b=b+1
          if(b.gt.NBF) then
           b=NA+1
           a=a+1
          endif
         enddo
         wmn(i,j,k)=wmn(i,j,k)*DSQRT(2.0d0)
         IF(TUNEMBPT) wmn2(i,j,k)=wmn2(i,j,k)*DSQRT(2.0d0) 
         if(i.ne.j) then
          wmn(j,i,k)=wmn(i,j,k)
          IF(TUNEMBPT) wmn2(j,i,k)=wmn2(i,j,k)
         endif
        enddo
       enddo 
      enddo
      end subroutine build_wmn

C  Loop to compute the EcGoWo energy. Recall that Go is used in GM EQ.
C  Very symple Eq: Ec = 2 sum _i sum_a sum_s [ (wia)^s ]**2 /(e(i)-e(a)-Omega(s))
      subroutine gw_gm_eq(NA,NCO,NBF,Nab,wmn,wmn2,EIG,EcGoWo,EcGMSOS,
     & TUNEMBPT,XpY,BIGOMEGA,ERImol)
      implicit none
      logical,intent(in)::TUNEMBPT
      integer,intent(in)::NCO,NBF,Nab,NA
      double precision,intent(inout)::EcGoWo,EcGMSOS
      double precision,dimension(Nab),intent(in)::BIGOMEGA
      double precision,dimension(NBF),intent(in)::EIG
      double precision,dimension(Nab,Nab),intent(in)::XpY
      double precision,dimension(NBF,NBF,Nab),intent(in)::wmn,wmn2
      double precision,dimension(NBF,NBF,NBF,NBF),intent(in)::ERImol
      integer::i,j,a,b,s,l,fst_virt
      double precision::tol10
      tol10=1.0D-10
      fst_virt=NA+1
      ! Doubly occ part (simplified version)
C      do a=NA+1,NBF
C       do i=1,NCO
C        do s=1,Nab
C         IF(TUNEMBPT) then
C         EcGoWo=EcGoWo+(wmn(i,a,s)*wmn2(i,a,s))/(EIG(i)-EIG(a)
C     &        -BIGOMEGA(s)+tol10)
C         ELSE
C         EcGoWo=EcGoWo+(wmn(i,a,s)**2.0d0)/(EIG(i)-EIG(a)
C     &        -BIGOMEGA(s)+tol10)
C         ENDIF
C        enddo
C       enddo
C      enddo
      ! Doubly occ part (4 index version)
       do a=NA+1,NBF
        do b=NA+1,NBF
         do i=1,NCO
          do j=1,NCO
           l=(j-1)*(NBF-NA)+(b-fst_virt+1)
           do s=1,Nab 
            EcGoWo=EcGoWo+wmn(i,a,s)*ERImol(j,i,a,b)*XpY(l,s)*
     &       DSQRT(2.0d0)/(EIG(i)-EIG(a)-BIGOMEGA(s)+tol10)
            EcGMSOS=EcGMSOS-wmn(i,a,s)*ERImol(j,i,b,a)*XpY(l,s)* ! Sign and ERI changed
     &       DSQRT(2.0d0)/(EIG(i)-EIG(a)-BIGOMEGA(s)+tol10)
           enddo
          enddo
         enddo
        enddo
       enddo
      ! Open-shell part
      if(NA/=NCO .and. TUNEMBPT) then  !TUNEMBPT=TRUE to have ERImol 
       do a=NA+1,NBF
        do b=NA+1,NBF
         do i=1,NCO
          do j=NCO+1,NA
           l=(j-1)*(NBF-NA)+(b-fst_virt+1)
           do s=1,Nab 
            EcGoWo=EcGoWo+0.5d0*wmn(i,a,s)*ERImol(j,i,a,b)*XpY(l,s)*
     &       DSQRT(2.0d0)/(EIG(i)-EIG(a)-BIGOMEGA(s)+tol10)
            EcGMSOS=EcGMSOS-0.5d0*wmn(i,a,s)*ERImol(j,i,b,a)*XpY(l,s)* ! Signed and ERI changed
     &       DSQRT(2.0d0)/(EIG(i)-EIG(a)-BIGOMEGA(s)+tol10)
           enddo
          enddo
         enddo
         do i=NCO+1,NA
          do j=1,NCO
           l=(j-1)*(NBF-NA)+(b-fst_virt+1)
           do s=1,Nab 
            EcGoWo=EcGoWo+0.5d0*wmn(i,a,s)*ERImol(j,i,a,b)*XpY(l,s)*
     &       DSQRT(2.0d0)/(EIG(i)-EIG(a)-BIGOMEGA(s)+tol10)
            EcGMSOS=EcGMSOS-0.5d0*wmn(i,a,s)*ERImol(j,i,b,a)*XpY(l,s)* ! Signed and ERI changed
     &       DSQRT(2.0d0)/(EIG(i)-EIG(a)-BIGOMEGA(s)+tol10)
           enddo
          enddo
          do j=NCO+1,NA
           if(i/=j) then
            l=(j-1)*(NBF-NA)+(b-fst_virt+1)
            do s=1,Nab 
             EcGoWo=EcGoWo+0.5d0*wmn(i,a,s)*ERImol(j,i,a,b)*XpY(l,s)*
     &        DSQRT(2.0d0)/(EIG(i)-EIG(a)-BIGOMEGA(s)+tol10)
             EcGMSOS=EcGMSOS-0.5d0*wmn(i,a,s)*ERImol(j,i,b,a)*XpY(l,s)* ! Signed and ERI changed
     &        DSQRT(2.0d0)/(EIG(i)-EIG(a)-BIGOMEGA(s)+tol10)
            enddo
           endif
          enddo
         enddo
        enddo
       enddo
      endif
      end subroutine gw_gm_eq

      subroutine mp2_eq(NA,NCO,NBF,EIG,ERImol,ERImol2,TUNEMBPT,EcMP2)
      implicit none
      logical,intent(in)::TUNEMBPT
      integer,intent(in)::NCO,NBF,NA
      double precision,intent(inout)::EcMP2
      double precision,dimension(NBF),intent(in)::EIG
      double precision,dimension(NBF,NBF,NBF,NBF),intent(in)::ERImol
      double precision,dimension(NBF,NBF,NBF,NBF),intent(in)::ERImol2
      integer::i,j,a,b
      double precision::tol10
      tol10=1.0D-10
      do a=NA+1,NBF
       do b=NA+1,NBF
        do i=1,NCO
         do j=1,NCO
          IF(TUNEMBPT) then
           ! Use [2Tij,ab - Tij,ba] as 2-RDM element (the tuned one in Piris PRA).
           EcMP2=EcMP2+ERImol2(a,b,j,i)*(2.0d0*ERImol(a,b,j,i)  
     &    -ERImol(a,b,i,j))/(EIG(i)+EIG(j)-EIG(a)-EIG(b)+tol10)
          ELSE
           EcMP2=EcMP2+ERImol(a,b,j,i)*(2.0d0*ERImol(a,b,j,i)
     &    -ERImol(a,b,i,j))/(EIG(i)+EIG(j)-EIG(a)-EIG(b)+tol10)
          ENDIF
         enddo
         do j=NCO+1,NA
          IF(TUNEMBPT) then
           ! Use [2Tij,ab - Tij,ba] as 2-RDM element (the tuned one in Piris PRA).
           EcMP2=EcMP2+ERImol2(a,b,j,i)*(ERImol(a,b,j,i)  
     &    -0.5D+0*ERImol(a,b,i,j))/(EIG(i)+EIG(j)-EIG(a)-EIG(b)+tol10)
          ELSE
           EcMP2=EcMP2+ERImol(a,b,j,i)*(ERImol(a,b,j,i)
     &    -0.5D+0*ERImol(a,b,i,j))/(EIG(i)+EIG(j)-EIG(a)-EIG(b)+tol10)
          ENDIF
         enddo
        enddo
        do i=NCO+1,NA
         do j=1,NCO
          IF(TUNEMBPT) then
           ! Use [2Tij,ab - Tij,ba] as 2-RDM element (the tuned one in Piris PRA).
           EcMP2=EcMP2+ERImol2(a,b,j,i)*(ERImol(a,b,j,i)  
     &    -0.5D+0*ERImol(a,b,i,j))/(EIG(i)+EIG(j)-EIG(a)-EIG(b)+tol10)
          ELSE
           EcMP2=EcMP2+ERImol(a,b,j,i)*(ERImol(a,b,j,i)
     &    -0.5D+0*ERImol(a,b,i,j))/(EIG(i)+EIG(j)-EIG(a)-EIG(b)+tol10)
          ENDIF
         enddo
         do j=NCO+1,NA
          IF(TUNEMBPT) then
           ! Use [2Tij,ab - Tij,ba] as 2-RDM element (the tuned one in Piris PRA).
           if(j/=i)EcMP2=EcMP2+0.5d+0*ERImol2(a,b,j,i)*(ERImol(a,b,j,i)
     &    -0.5D+0*ERImol(a,b,i,j))/(EIG(i)+EIG(j)-EIG(a)-EIG(b)+tol10)
          ELSE
           if(j/=i)EcMP2=EcMP2+0.5d+0*ERImol(a,b,j,i)*(ERImol(a,b,j,i)
     &    -0.5D+0*ERImol(a,b,i,j))/(EIG(i)+EIG(j)-EIG(a)-EIG(b)+tol10)
          ENDIF
         enddo
        enddo
       enddo
      enddo 
      end subroutine mp2_eq

      subroutine rpa_sosex_eq(NA,NCO,NBF,Nab,ERImol,ERImol2,EIG,
     & BIGOMEGA,wmn,iEcRPA,iEcSOSEX,TUNEMBPT)
      implicit none
      logical,intent(in)::TUNEMBPT
      integer,intent(in)::NCO,NBF,Nab,NA
      double precision,intent(inout)::iEcRPA,iEcSOSEX
      double precision,dimension(Nab),intent(in)::BIGOMEGA
      double precision,dimension(NBF),intent(in)::EIG
      double precision,dimension(NBF,NBF,Nab),intent(in)::wmn
      double precision,dimension(NBF,NBF,NBF,NBF),intent(in)::ERImol
      double precision,dimension(NBF,NBF,NBF,NBF),intent(in)::ERImol2
      double precision,allocatable,dimension(:)::freqs,weights
      complex,allocatable,dimension(:)::cfreqs
      integer::order,i,j,a,b
      double precision::PI,r0,rlast,integral,integrated_omega
      PI=3.141592653589793238D+00
C  Quadrature 0 to Inf
      r0=0.0d0
      rlast=1.0d0
      order=40 ! we can keep it fix it, FHI-aims guys were right this time!!!
      ALLOCATE(weights(order),freqs(order),cfreqs(order))
      call gauss_legendre(r0,rlast,order)
      open(unit=314,file='tempq_x.txt',status='old')
      open(unit=315,file='tempq_w.txt',status='old')
      do i=1,order
       read(314,*) freqs(i)
       read(315,*) weights(i)
       weights(i)=weights(i)/(1.0d0-freqs(i))**2.0d0
       freqs(i)=freqs(i)/(1.0d0-freqs(i))
       cfreqs(i)=cmplx(0.0d0,freqs(i))
      enddo 
      close(314)
      close(315)
      call system('/bin/rm -r tempq*txt')
C  RPA + SOSEX integrated
      do a=NA+1,NBF
       do b=NA+1,NBF
        do i=1,NCO
         do j=1,NCO
          integral=integrated_omega(i,j,a,b,order,NBF,Nab,wmn,weights,
     &    freqs,cfreqs,EIG,BIGOMEGA,ERImol(i,j,b,a))
          IF(TUNEMBPT) THEN
           iEcRPA=iEcRPA-ERImol2(i,j,b,a)*integral      ! Notice the minus, stupid FHI-aims developers.
           iEcSOSEX=iEcSOSEX+ERImol2(i,j,a,b)*integral  ! Opposite sign of RPA and dif. 2e- integral
          ELSE
           iEcRPA=iEcRPA-ERImol(i,j,b,a)*integral       ! Notice the minus, stupid FHI-aims developers.
           iEcSOSEX=iEcSOSEX+ERImol(i,j,a,b)*integral   ! Opposite sign of RPA and dif. 2e- integral
          ENDIF
         enddo
         do j=NCO+1,NA
          integral=integrated_omega(i,j,a,b,order,NBF,Nab,wmn,weights,
     &    freqs,cfreqs,EIG,BIGOMEGA,ERImol(i,j,b,a))
          IF(TUNEMBPT) THEN
           iEcRPA=iEcRPA-0.5d+0*ERImol2(i,j,b,a)*integral      ! Notice the minus, stupid FHI-aims developers.
           iEcSOSEX=iEcSOSEX+0.5d+0*ERImol2(i,j,a,b)*integral  ! Opposite sign of RPA and dif. 2e- integral
          ELSE
           iEcRPA=iEcRPA-0.5d+0*ERImol(i,j,b,a)*integral       ! Notice the minus, stupid FHI-aims developers.
           iEcSOSEX=iEcSOSEX+0.5d+0*ERImol(i,j,a,b)*integral   ! Opposite sign of RPA and dif. 2e- integral
          ENDIF
         enddo
        enddo
        do i=NCO+1,NA
         do j=1,NCO
          integral=integrated_omega(i,j,a,b,order,NBF,Nab,wmn,weights,
     &    freqs,cfreqs,EIG,BIGOMEGA,ERImol(i,j,b,a))
          IF(TUNEMBPT) THEN
           iEcRPA=iEcRPA-0.5d+0*ERImol2(i,j,b,a)*integral      ! Notice the minus, stupid FHI-aims developers.
           iEcSOSEX=iEcSOSEX+0.5d+0*ERImol2(i,j,a,b)*integral  ! Opposite sign of RPA and dif. 2e- integral
          ELSE
           iEcRPA=iEcRPA-0.5d+0*ERImol(i,j,b,a)*integral       ! Notice the minus, stupid FHI-aims developers.
           iEcSOSEX=iEcSOSEX+0.5d+0*ERImol(i,j,a,b)*integral   ! Opposite sign of RPA and dif. 2e- integral
          ENDIF
         enddo
         do j=NCO+1,NA
          integral=integrated_omega(i,j,a,b,order,NBF,Nab,wmn,weights,
     &    freqs,cfreqs,EIG,BIGOMEGA,ERImol(i,j,b,a))
          IF(TUNEMBPT) THEN
           if(j/=i)iEcRPA=iEcRPA-0.5d+0*ERImol2(i,j,b,a)*integral      ! Notice the minus, stupid FHI-aims developers.
           if(j/=i)iEcSOSEX=iEcSOSEX+0.5d+0*ERImol2(i,j,a,b)*integral  ! Opposite sign of RPA and dif. 2e- integral
          ELSE
           if(j/=i)iEcRPA=iEcRPA-0.5d+0*ERImol(i,j,b,a)*integral       ! Notice the minus, stupid FHI-aims developers.
           if(j/=i)iEcSOSEX=iEcSOSEX+0.5d+0*ERImol(i,j,a,b)*integral   ! Opposite sign of RPA and dif. 2e- integral
          ENDIF
         enddo
        enddo
       enddo
      enddo
      iEcRPA=iEcRPA/PI ! Assuming int 0 to 1 W_lambda(iw) d_lambda approx W_1(iw)/2 (Trapezoid rule!), also x 4 (4 is for spins of X_o, or 2e- integrals)
      iEcSOSEX=0.5d+0*iECSOSEX/PI ! Use 1/2 because spinless 2e- integrals (only alpha,alpha,alpha,alpha and beta,beta,beta,beta)
      DEALLOCATE(weights,freqs,cfreqs) 
      end subroutine rpa_sosex_eq
      
      subroutine ccsd_eq(NA,NCO,NBF,ERImol,ERImol2,EIG,EcCCSD,TUNEMBPT,
     & CCSD_READ)
      USE m_ccsd
      implicit none
      logical,intent(in)::TUNEMBPT,CCSD_READ
      integer,intent(in)::NCO,NBF,NA
      double precision,intent(inout)::EcCCSD
      double precision,dimension(NBF),intent(in)::EIG
      double precision,dimension(NBF,NBF,NBF,NBF),intent(inout)::ERImol
      double precision,dimension(NBF,NBF,NBF,NBF),intent(inout)::ERImol2
      double precision,allocatable,dimension(:,:)::FockM
      double precision,allocatable,dimension(:,:,:,:)::ERItmp
      integer::p,q,r,s,Nocc,verbose=0,iter=0,NCO2,NA2,maxiter=1000000
      double precision::error=1.0d0,Eccsd,Eccsd_new,tol7=1.0d-7
      allocate(FockM(NBF,NBF),ERItmp(NBF,NBF,NBF,NBF))
      FockM=0.0d0
      do p=1,NBF
       FockM(p,p)=EIG(p)     
      enddo
      ERItmp=0.0d0
      do p=1,NBF
       do q=1,NBF
        do r=1,NBF
         do s=1,NBF
          if(abs(ERImol(p,q,s,r))>tol7) then
           ERItmp(p,q,r,s)=ERImol(p,q,s,r)
          endif
         enddo
        enddo
       enddo
      enddo
      ERImol=0.0d0
      ERImol=ERItmp
      IF(TUNEMBPT) then
       ERItmp=0.0d0
       do p=1,NBF
        do q=1,NBF
         do r=1,NBF
          do s=1,NBF
           if(abs(ERImol2(p,q,s,r))>tol7) then
            ERItmp(p,q,r,s)=ERImol2(p,q,s,r)
           endif
          enddo
         enddo
        enddo
       enddo
       ERImol2=0.0d0
       ERImol2=ERItmp
      ENDIF
      deallocate(ERItmp)
      if(NA-NCO==0) then
       Nocc=NCO
       NCO2=NCO*2
       NA2=NCO2
      else
       Nocc=NCO+NA
       NCO2=NCO*2
       NA2=NA*2
      endif
      call ccsd_init(NBF,Nocc,verbose,FockM,ERImol) 
      if(CCSD_READ) then
       call ccsd_read_guess() 
       write(*,*) ' '
       write(*,*) ' Read the T amplitudes file'
       write(*,*) ' '
      endif
      deallocate(FockM)
      write(*,*) ' '
      do
       if(error>tol7 .and. iter<maxiter) then
        call ccsd_update_ts(ERImol)
        IF(TUNEMBPT) then
         Eccsd_new=ccsd_en_nof(NCO2,NA2,ERImol2)
        ELSE
         Eccsd_new=ccsd_en_nof(NCO2,NA2,ERImol)
        ENDIF
        error=abs(Eccsd-Eccsd_new)
        Eccsd=Eccsd_new
        iter=iter+1
        if(mod(iter,2)==0) then
         write(*,'(a,F15.10,a,i5)') ' CCSD energy ',Eccsd,' iter. ',iter
        endif
       else
        exit
       endif
      enddo
      write(*,*) ' '
      write(*,*) 'CCSD iterative procedure finished after ',iter,
     & ' iterations'
      write(*,*) ' '
      EcCCSD=Eccsd_new
      call ccsd_write_last()
      call ccsd_clean()
      end subroutine ccsd_eq

      subroutine transform2(NBF,TEMPM2,ERImol,ERImol2)
      implicit none
      integer,intent(in)::NBF
      double precision,dimension(nbf,nbf,nbf,nbf),intent(inout)::ERImol
      double precision,dimension(nbf,nbf,nbf,nbf),intent(inout)::ERImol2
      double precision,dimension(nbf,nbf),intent(in)::TEMPM2
      integer::i,j,k,l,m
      double precision::tol10,aux,aux2
      tol10=1.0D-10
      write(*,*) ' '
      write(*,*) ' Starting transformation <IJ|KL> -> <PQ|RS>'
      ! L -> S
      open(unit=66,form='unformatted')
      open(unit=67,form='unformatted')
      do i=1,NBF
       do j=1,NBF
        do k=1,NBF
         do m=1,NBF
          aux=0.0d0
          aux2=0.0d0
          do l=1,NBF
           aux=aux+TEMPM2(l,m)*ERImol(i,j,k,l)
           aux2=aux2+TEMPM2(l,m)*ERImol2(i,j,k,l)
          enddo
          if(abs(aux).gt.tol10)  write(66) i,j,k,m,aux
          if(abs(aux2).gt.tol10) write(67) i,j,k,m,aux2
         enddo
        enddo
       enddo
      enddo
      close(66)
      close(67)
      ERImol=0.0d0
      ERImol2=0.0d0
      open(unit=66,form='unformatted',status='old')
       do while(.true.)
        read(66,end=99) i,j,k,l,aux
        if (i.eq.0) goto 99
         ERImol(i,j,k,l)=aux
       enddo
   99  continue
      close(66)
      open(unit=67,form='unformatted',status='old')
       do while(.true.)
        read(67,end=100) i,j,k,l,aux
        if (i.eq.0) goto 100
         ERImol2(i,j,k,l)=aux
       enddo
  100 continue
      close(67)
      call system("/bin/rm fort.66 fort.67")
      write(*,*) ' L -> S done'
      ! K -> R
      open(unit=66,form='unformatted')
      open(unit=67,form='unformatted')
      do i=1,NBF
       do j=1,NBF
        do m=1,NBF
         do l=1,NBF
          aux=0.0d0
          aux2=0.0d0
          do k=1,NBF
           aux=aux+TEMPM2(k,m)*ERImol(i,j,k,l)
           aux2=aux2+TEMPM2(k,m)*ERImol2(i,j,k,l)
          enddo
          if(abs(aux).gt.tol10)  write(66) i,j,m,l,aux
          if(abs(aux2).gt.tol10) write(67) i,j,m,l,aux2
         enddo
        enddo
       enddo
      enddo
      close(66)
      close(67)
      ERImol=0.0d0
      ERImol2=0.0d0
      open(unit=66,form='unformatted',status='old')
       do while(.true.)
        read(66,end=101) i,j,k,l,aux
        if (i.eq.0) goto 101
         ERImol(i,j,k,l)=aux
       enddo
  101 continue
      close(66)
      open(unit=67,form='unformatted',status='old')
       do while(.true.)
        read(67,end=102) i,j,k,l,aux
        if (i.eq.0) goto 102
         ERImol2(i,j,k,l)=aux
       enddo
  102 continue
      close(67)
      call system("/bin/rm fort.66 fort.67")
      write(*,*) ' K -> R done'
      ! J -> Q
      open(unit=66,form='unformatted')
      open(unit=67,form='unformatted')
      do i=1,NBF
       do m=1,NBF
        do k=1,NBF
         do l=1,NBF
          aux=0.0d0
          aux2=0.0d0
          do j=1,NBF
           aux=aux+TEMPM2(j,m)*ERImol(i,j,k,l)
           aux2=aux2+TEMPM2(j,m)*ERImol2(i,j,k,l)
          enddo
          if(abs(aux).gt.tol10)  write(66) i,m,k,l,aux
          if(abs(aux2).gt.tol10) write(67) i,m,k,l,aux2
         enddo
        enddo
       enddo
      enddo
      close(66)
      close(67)
      ERImol=0.0d0
      ERImol2=0.0d0
      open(unit=66,form='unformatted',status='old')
       do while(.true.)
        read(66,end=103) i,j,k,l,aux
        if (i.eq.0) goto 103
         ERImol(i,j,k,l)=aux
       enddo
  103 continue
      close(66)
      open(unit=67,form='unformatted',status='old')
       do while(.true.)
        read(67,end=104) i,j,k,l,aux
        if (i.eq.0) goto 104
         ERImol2(i,j,k,l)=aux
       enddo
  104 continue
      close(67)
      call system("/bin/rm fort.66 fort.67")
      write(*,*) ' J -> Q done'
      ! I -> P
      open(unit=66,form='unformatted')
      open(unit=67,form='unformatted')
      do m=1,NBF
       do j=1,NBF
        do k=1,NBF
         do l=1,NBF
          aux=0.0d0
          aux2=0.0d0
          do i=1,NBF
           aux=aux+TEMPM2(i,m)*ERImol(i,j,k,l)
           aux2=aux2+TEMPM2(i,m)*ERImol2(i,j,k,l)
          enddo
          if(abs(aux).gt.tol10)  write(66) m,j,k,l,aux
          if(abs(aux2).gt.tol10) write(67) m,j,k,l,aux2
         enddo
        enddo
       enddo
      enddo
      close(66)
      close(67)
      ERImol=0.0d0
      ERImol2=0.0d0
      open(unit=66,form='unformatted',status='old')
       do while(.true.)
        read(66,end=105) i,j,k,l,aux
        if (i.eq.0) goto 105
         ERImol(i,j,k,l)=aux
       enddo
  105 continue
      close(66)
      open(unit=67,form='unformatted',status='old')
       do while(.true.)
        read(67,end=106) i,j,k,l,aux
        if (i.eq.0) goto 106
         ERImol2(i,j,k,l)=aux
       enddo
  106 continue
      close(67)
      call system("/bin/rm fort.66 fort.67")
      write(*,*) ' I -> P done'
      write(*,*) ' Transformation finished'
      write(*,*) ' '
      end subroutine transform2

      subroutine transform2mem(NBF,TEMPM2,ERImol)
      implicit none
      integer,intent(in)::NBF
      double precision,dimension(nbf,nbf,nbf,nbf),intent(inout)::ERImol
      double precision,dimension(nbf,nbf),intent(in)::TEMPM2
      double precision,dimension(:,:,:,:),allocatable::TMPDM2
      integer::i,j,k,l,m
      write(*,*) ' '
      write(*,*) ' Starting transformation <IJ|KL> -> <PQ|RS>'
      ! L -> S
      allocate(TMPDM2(NBF,NBF,NBF,NBF))
      do i=1,NBF
       do j=1,NBF
        do k=1,NBF
         do m=1,NBF
          TMPDM2(i,j,k,m)=0.0d0
          do l=1,NBF
           TMPDM2(i,j,k,m)=TMPDM2(i,j,k,m)+TEMPM2(l,m)*ERImol(i,j,k,l)
          enddo
         enddo
        enddo
       enddo
      enddo
      write(*,*) ' L -> S done'
      ! K -> R
      do i=1,NBF
       do j=1,NBF
        do m=1,NBF
         do l=1,NBF
          ERImol(i,j,m,l)=0.0d0
          do k=1,NBF
           ERImol(i,j,m,l)=ERImol(i,j,m,l)+TEMPM2(k,m)*TMPDM2(i,j,k,l)
          enddo
         enddo
        enddo
       enddo
      enddo
      write(*,*) ' K -> R done'
      ! J -> Q
      do i=1,NBF
       do m=1,NBF
        do k=1,NBF
         do l=1,NBF
          TMPDM2(i,m,k,l)=0.0d0
          do j=1,NBF
           TMPDM2(i,m,k,l)=TMPDM2(i,m,k,l)+TEMPM2(j,m)*ERImol(i,j,k,l)
          enddo
         enddo
        enddo
       enddo
      enddo
      write(*,*) ' J -> Q done'
      ! I -> P
      do m=1,NBF
       do j=1,NBF
        do k=1,NBF
         do l=1,NBF
          ERImol(m,j,k,l)=0.0d0
          do i=1,NBF
           ERImol(m,j,k,l)=ERImol(m,j,k,l)+TEMPM2(i,m)*TMPDM2(i,j,k,l)
          enddo
         enddo
        enddo
       enddo
      enddo
      write(*,*) ' I -> P done'
      write(*,*) ' Transformation finished'
      write(*,*) ' '
      deallocate(TMPDM2)
      end subroutine transform2mem

C      subroutine print2eint(nbf,erimol)
C      implicit none
C      integer,intent(in)::nbf
C      double precision,dimension(nbf,nbf,nbf,nbf),intent(in)::erimol
C      integer::i,j,k,l
C      do i=1,nbf
C       do j=1,nbf
C        do k=1,nbf
C         do l=1,nbf 
C           write(*,*) i,j,k,l,ERImol(i,j,k,l)
C         enddo
C        enddo
C       enddo
C      enddo
C      end subroutine print2eint

!======================================================================!
