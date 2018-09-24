 ! THE BACKGROUND SOLVER. 

    MODULE BACKGROUND_PATH
       USE RKSUITE_90_PREC

       REAL(WP),DIMENSION(:),ALLOCATABLE::EF,PHI,PHID,PHIDD,HUBBLE

       INTEGER::ESEARCH

       REAL(WP)::AFI,AI,NEND,H_INI,END_OF_INFLATION,MPE
       
       Logical::force_ai,multiphase
             

     END MODULE BACKGROUND_PATH



     SUBROUTINE BACKGROUND_VALUE
       USE RKSUITE_90_PREC
       USE RKSUITE_90
       USE THEORY_PARAM
       USE PARAMETERS
       USE BACKGROUND_PATH
       USE INTERFACING, ONLY : POTENTIAL_POTENTIALPRIME,BACKGROUND
       IMPLICIT NONE



       REAL(WP),ALLOCATABLE,DIMENSION(:)::BGND_INI,POTENTIAL,POTENTIAL_PRIME,XX,Y_GOT,YDERIV_GOT,THRESHOLDS,yprime
       REAL(WP),ALLOCATABLE,DIMENSION(:,:)::Y
       REAL(WP)::STEPSIZE,E_H
       REAL(WP)::POT_INI,PHID_INI,TOLERANCE
       INTEGER::STEPS,I,II,JJ,FLAG
       REAL(WP)::NEND_50,HUBBLE_AT_NEND_50,POT_PHI_INI,NSTART,NSTOP,N_GOT
       INTEGER,PARAMETER::FILE_NUM=20,FILE_NUM1=21,FILE_NUM2=22,FILE_NUM3=23
       TYPE(RK_COMM_REAL_1D) :: COMMB

       !===================DECLARATION=======================



       STEPSIZE=(NF-NI)/(DBLE(BGSTEPS))

       
       ALLOCATE(BGND_INI(2))


       !=====================================================
        ALLOCATE(XX(BGSTEPS),Y(2,BGSTEPS),Y_GOT(2),YDERIV_GOT(2),YPRIME(BGSTEPS),THRESHOLDS(2))
        TOLERANCE=1.0D-6
        THRESHOLDS=1.0D-8
       !===============BACKGROUND INITIAL CONDITION=========

       BGND_INI(1)=INITIALPHI!INITIAL PHI


       !-------------THE POTENTIAL DEFINED HERE TOO---------
       !TO GET THE HUBBLE PARAMETER INITIAL VALUE

       CALL POTENTIAL_POTENTIALPRIME(POT_INI,POT_PHI_INI,BGND_INI(1))


       !----------------------------------------------------

       PHID_INI=INITIALPHIDOT
       H_INI=(1.0_WP/3.0_WP)*((0.5_WP*(PHID_INI**2.0_WP))+POT_INI)
       H_INI=dSQRT(H_INI)
       BGND_INI(2)=PHID_INI/H_INI

       PRINT*,'THE INITIAL CONDITION FOR THE INFLATIONS ARE'
       PRINT*,'\PHI_I=',BGND_INI(1)
       PRINT*,'\PHI_NI',BGND_INI(2)


  
  !================== STORING ZEROTH VALUES ======================================
    XX(1)=NI
    Y(:,1)=BGND_INI    
    yprime(1)=-(3.0_wp-(0.5_wp*BGND_INI(2)*BGND_INI(2)))*BGND_INI(2)-POT_PHI_INI/(H_INI*H_INI)
  !===============================================================================
   
#ifndef MPI
  OPEN(UNIT=FILE_NUM,FILE='plots/EPSILON.dat',FORM='FORMATTED',STATUS='REPLACE')
#endif
  CALL SETUP(COMMB,NI,BGND_INI,NF,TOLERANCE,THRESHOLDS,method="M",MESSAGE=.TRUE.)

   ADAPTIVE_BGRND:DO I=1,(BGSTEPS-1)

     NSTART=NI+((DBLE(I)-1.0_WP)*(NF-NI)/DBLE(BGSTEPS-1.0_WP))
     NSTOP=NI+(DBLE(I)*(NF-NI)/DBLE(BGSTEPS-1.0_WP))



     CALL range_integrate(COMMB,BACKGROUND,NSTOP,N_GOT,Y_GOT,YDERIV_GOT,FLAG=FLAG)
          IF(FLAG/=1) THEN
            PRINT*,"INTEGRATION NOT SUCCESSFUL"
            PRINT*,'COULD EVALUATE TILL N=',N_GOT,"INSTEAD N=",NSTOP
            IF(ABS(N_GOT-Nstop)>1.0d-10) STOP "EXITING"
          ENDIF
     !==================== STORE INTERMEDIATE VALUES TILL NEND ==================
     XX(I+1)=NSTOP
     Y(:,I+1)=Y_GOT(:)
     yprime(I+1)=YDERIV_GOT(2)
     STEPS=I+1
     !===========================================================================
     
     !DEFINING THE EPSILON
     
     E_H=0.0_WP!SETTING \EPSILON_H=0 IN EACH STEP

     E_H=0.5_WP*Y(2,I+1)*Y(2,I+1)
     !CHECK FOR THE END OF INFLATION 

     IF(.not.(multiphase).and.E_H>1.0_WP) THEN

      END_OF_INFLATION=XX(I+1)
      elseif (multiphase.and.E_H>1.0_WP.and.XX(I+1)>MPE) then
      END_OF_INFLATION=XX(I+1)

     ENDIF
     
     IF(.not.(multiphase).and.E_H>1.0_WP) EXIT
     if(multiphase.and.E_H>1.0_WP.and.XX(I+1)>MPE) EXIT
#ifndef MPI    
     WRITE(FILE_NUM,*),XX(I+1),E_H
#endif
     IF(I==BGSTEPS-1.and.E_H<1.0_WP) THEN        
       END_OF_INFLATION=XX(I+1)
       PRINT*,"CAUTION : INFLATION HAS NOT ENDED: EXITING AT FINAL EFOLD PROVIDED"
     ENDIF  
     
     IF(I==BGSTEPS-1.and.E_H<1.0_WP)  EXIT 

  END DO ADAPTIVE_BGRND
#ifndef MPI
  CLOSE(FILE_NUM)
#endif
  CALL COLLECT_GARBAGE(COMMB)
 
  !======== STEPS TO THE END OF INFLATION ========================================
  
  
       ESEARCH=STEPS


       !==========================================================
       IF(ALLOCATED(HUBBLE))DEALLOCATE(HUBBLE)
       IF(ALLOCATED(EF))DEALLOCATE(EF)
       IF(ALLOCATED(PHI))DEALLOCATE(PHI)
       IF(ALLOCATED(PHID))DEALLOCATE(PHID)
       IF(ALLOCATED(PHIDD))DEALLOCATE(PHIDD)

       ALLOCATE(HUBBLE(STEPS),POTENTIAL(STEPS),POTENTIAL_PRIME(STEPS),EF(STEPS),PHI(STEPS),PHID(STEPS),PHIDD(STEPS))

       EF(1:STEPS)=XX(1:STEPS)
  
       PHI(1:STEPS)=Y(1,1:STEPS)

       PHID(1:STEPS)=Y(2,1:STEPS)

       PHIDD(1:STEPS)=yprime(1:STEPS)
       !------------------------------------------------------
       STORE_LOOP:DO II=1,STEPS
        
          CALL POTENTIAL_POTENTIALPRIME(POTENTIAL(II),POTENTIAL_PRIME(II),PHI(II))

          HUBBLE(II)=dSQRT(POTENTIAL(II)/(3.0_WP-0.5_WP*(PHID(II)**2.0_WP)))


       END DO STORE_LOOP


! 
       !================== DEALLOCATE THE USELESS MEMORY ====================
      DEALLOCATE(XX,Y,Y_GOT,YDERIV_GOT,THRESHOLDS,yprime)
       !-----------------------------------------------------
#ifndef MPI

     OPEN(UNIT=FILE_NUM,FILE='plots/PHI.dat',FORM='FORMATTED',STATUS='REPLACE')
     OPEN(UNIT=FILE_NUM1,FILE='plots/DPHI_DN.dat',FORM='FORMATTED',STATUS='REPLACE')
     OPEN(UNIT=FILE_NUM2,FILE='plots/HUBBLE.dat',FORM='FORMATTED',STATUS='REPLACE')
     OPEN(UNIT=FILE_NUM3,FILE='plots/PHINN.dat',FORM='FORMATTED',STATUS='REPLACE')
     
     KEEP_DATA:DO II=1,STEPS
        WRITE(FILE_NUM,*)EF(II),PHI(II)
        WRITE(FILE_NUM1,*)EF(II),PHID(II)
        WRITE(FILE_NUM2,*)EF(II),HUBBLE(II)
        WRITE(FILE_NUM3,*)EF(II),PHIDD(II)
     END DO KEEP_DATA
     

  
     CLOSE(FILE_NUM)
     CLOSE(FILE_NUM1) 
     CLOSE(FILE_NUM2)
     CLOSE(FILE_NUM3)
#endif
        NEND_50=END_OF_INFLATION-N_pivot
        NEND=EF(STEPS)
        
        PRINT*,'INFLATION ENDS AT EFOLDS=',END_OF_INFLATION

        NEND_50_LOOP:DO JJ=1,STEPS
           HUBBLE_AT_NEND_50=HUBBLE(JJ)
           IF(EF(JJ)>NEND_50)EXIT
        END DO NEND_50_LOOP

        IF(FORCE_AI) THEN 
          AI=AFI
          PRINT*,'IMPOSED A_I=',AI
        ELSE
        !-------------COMPUTING AI-------------------------------------------
          AI=PIVOTK/(HUBBLE_AT_NEND_50*DEXP(NEND_50))
          PRINT*,'CALCULATED A_I=',AI
       ENDIF
       
       PRINT*,'============================================================================='
       DEALLOCATE(POTENTIAL,POTENTIAL_PRIME)

     END SUBROUTINE BACKGROUND_VALUE
