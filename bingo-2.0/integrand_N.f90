!  THE INTEGRAND FOR G_4. CHANGE THE factor AND THE integrand TO CALCULATE OTHER G  

   FUNCTION INTEGRAND(NFOLD)
      
       USE THEORY_PARAM
       USE PARAMETERS
       USE BACKGROUND_PATH
       USE ODE_K
       USE POT2DER
       USE STOREK
       USE INITIATION
       USE MODESTORE
       USE FNLTERM
       USE RKSUITE_90_PREC
       USE RKSUITE_90
       USE INTERFACING, ONLY : POTENTIAL_POTENTIALPRIME,INITIAL_COND,CURV_PERT,FIND_NIC
       IMPLICIT NONE
       

       REAL(WP),INTENT(IN)::NFOLD
       REAL(WP),ALLOCATABLE,DIMENSION(:)::PERT_INI,Y_GOT,YDERIV_GOT,THRESHOLDS
       REAL(WP)::TOLERANCE,NFOLD_GOT,cutoff
       REAL(WP)::HUBB_N,EPSILON_ETAPRIME,EPSILON1
       REAL(WP)::HUBB,INTEGRAND,PHI_N,PHI_N_N,PHI_N_N_N,VTOT,V_PHI
       COMPLEX(WP)::CURV_K1,CURV_PRIME_K1,CURV_K2,CURV_PRIME_K2,CURV_K3,CURV_PRIME_K3,FACTOR
       COMPLEX(WP)::INTEGRAND1,INTEGRAND2,INTEGRAND3,INTEGRAND4,INTEGRAND5,INTEGRAND6
       REAL(WP),PARAMETER::TINYVALUE=0.1d0
       TYPE(RK_COMM_REAL_1D) :: COMM
       INTEGER :: FLAG

       IF((NFOLD<NIC).OR.(NFOLD>NEND)) STOP 'INTEGRATION LIMIT IS LESS THAN THE HUBBLE EXIT '
        
  
       !===================ALLOCATION========================
       ALLOCATE(PERT_INI(6),Y_GOT(6),YDERIV_GOT(6),THRESHOLDS(6))

       TOLERANCE=(1D-5)*(10.0**(-acc))
       THRESHOLDS(1:6)=(1.0D-7)*(10.0**(-acc))
       !================= INITIATION OF THE INTEGRAND =======

       INTEGRAND=0.0_WP
       PHI_N_N=0.0_WP
       PHI_N_N_N=0.0_WP
       FACTOR=0.0_WP
       HUBB=0.0_WP
       
       CALL FIND_NIC(modek1)
       
       CALL INITIAL_COND(MODEK1)
       
       MODEK=MODEK1

       
       !============QUANTITY AT SUB HUBBLE SCALES===========================
       PERT_INI=0.0_WP
       PERT_INI(1)=PHI_IC
       PERT_INI(2)=PHID_IC
       PERT_INI(3)=R_INI !REAL PART OF R_K
       PERT_INI(4)=RPRIME_RE !REAL PART OF R_K PRIME
       PERT_INI(5)=0.0_WP!IMAGIMARY PART OF R_K
       PERT_INI(6)=RPRIME_IM!IMAGINARY PART OF R_K PRIME
          
          

       
       IF(NFOLD>(NIC+1.0D-10)) THEN

          CALL SETUP(COMM,NIC,PERT_INI,END_OF_INFLATION,TOLERANCE,THRESHOLDS,Method="H",MESSAGE=.FALSE.)
          CALL RANGE_INTEGRATE(COMM,CURV_PERT,NFOLD,NFOLD_GOT,Y_GOT,YDERIV_GOT,FLAG=FLAG)
          IF(FLAG/=1) THEN
            PRINT*,"INTEGRATION NOT SUCCESSFUL"
            PRINT*,'COULD EVALUATE TILL N=',NFOLD_GOT,"INSTEAD N=",NFOLD
            IF(ABS(NFOLD_GOT-NFOLD)>1.0d-10) STOP "EXITING"
          ENDIF
          PERT_INI(:)=Y_GOT(:)
       ENDIF


        CALL POTENTIAL_POTENTIALPRIME(VTOT,V_PHI,PERT_INI(1))
  
       
       HUBB=dSQRT(VTOT/(3.0_WP-0.5_WP*(PERT_INI(2)**2.0_WP)))

       PHI_N=PERT_INI(2)

      

       !=============   DEFINTION OF RK1 AND RK1' ==================================== 


       CURV_K1=CMPLX(PERT_INI(3),PERT_INI(5),KIND=WP)

       
       CURV_PRIME_K1=CMPLX(PERT_INI(4),PERT_INI(6),KIND=WP)

       !============================================================================
     
       HUBB_N=(-0.5D0*HUBB*PHI_N*PHI_N)

       PHI_N_N=-V_PHI/(HUBB*HUBB)-(3.0D0-(0.5D0*PHI_N*PHI_N))*PHI_N

       PHI_N_N_N=(-V_PHI_PHI*PHI_N/(HUBB*HUBB))+(2.0D0*V_PHI*HUBB_N/(HUBB*HUBB*HUBB))
       PHI_N_N_N=PHI_N_N_N-(3.0D0*PHI_N_N)+1.5D0*PHI_N*PHI_N*PHI_N_N      

       EPSILON1=0.5D0*PHI_N*PHI_N
       EPSILON_ETAPRIME=AI*dEXP(NFOLD)*HUBB*((PHI_N_N_N*PHI_N)-PHI_N_N*PHI_N_N)
       
       EPS2=2.0D0*PHI_N_N/PHI_N
       
      CALL COLLECT_GARBAGE(COMM)      
!=======================  changed =======================================
      ! If equilateral
     if(equilateral) then 
	CURV_K2=CURV_K1
	CURV_K3=CURV_K1
	CURV_PRIME_K2=CURV_PRIME_K1
	CURV_PRIME_K3=CURV_PRIME_K1 
     end if

      !If Isosceles or scalene
     integration_for_isosceles_or_scalene:if(isosceles.or.scalene) then
     if(isosceles) then
	CURV_K2=CURV_K1
	CURV_PRIME_K2=CURV_PRIME_K1
     end if 
        CALL FIND_NIC(modek3)
        call initial_cond(modek3)
        MODEK=MODEK3
     !============QUANTITY AT SUB HUBBLE SCALES===========================
     !========= These initial conditions shall be changed according to modek3 ============
       PERT_INI=0.0_WP
       PERT_INI(1)=PHI_IC
       PERT_INI(2)=PHID_IC
       PERT_INI(3)=R_INI !REAL PART OF R_K
       PERT_INI(4)=RPRIME_RE !REAL PART OF R_K PRIME
       PERT_INI(5)=0.0_WP!IMAGIMARY PART OF R_K
       PERT_INI(6)=RPRIME_IM!IMAGINARY PART OF R_K PRIME
     !===================================================================== 
       IF(NFOLD>(NIC+1.0D-10)) THEN

          CALL SETUP(COMM,NIC,PERT_INI,END_OF_INFLATION,TOLERANCE,THRESHOLDS,Method="H",MESSAGE=.FALSE.)
          CALL RANGE_INTEGRATE(COMM,CURV_PERT,NFOLD,NFOLD_GOT,Y_GOT,YDERIV_GOT,FLAG=FLAG)
          IF(FLAG/=1) THEN
            PRINT*,"INTEGRATION NOT SUCCESSFUL"
            PRINT*,'COULD EVALUATE TILL N=',NFOLD_GOT,"INSTEAD N=",NFOLD,"In modek3=",modek3
            IF(ABS(NFOLD_GOT-NFOLD)>1.0d-10) STOP "EXITING"
          ENDIF
          PERT_INI(:)=Y_GOT(:)
       ENDIF
       
       
       !=============   DEFINTION OF RK3 AND RK3' ==================================== 


       CURV_K3=CMPLX(PERT_INI(3),PERT_INI(5),KIND=WP)

       
       CURV_PRIME_K3=CMPLX(PERT_INI(4),PERT_INI(6),KIND=WP)

       
       CALL COLLECT_GARBAGE(COMM)      

      !=========  At nfold background quantities shall not chage due to change of mode ========
    end if integration_for_isosceles_or_scalene

      !If scalene only

     integration_for_scalene:if(scalene) then
        CALL FIND_NIC(modek2)
        call initial_cond(modek2)
        MODEK=MODEK2
     !============QUANTITY AT SUB HUBBLE SCALES===========================
     !========= These initial conditions shall be changed according to modek2 ============
       PERT_INI=0.0_WP
       PERT_INI(1)=PHI_IC
       PERT_INI(2)=PHID_IC
       PERT_INI(3)=R_INI !REAL PART OF R_K
       PERT_INI(4)=RPRIME_RE !REAL PART OF R_K PRIME
       PERT_INI(5)=0.0_WP!IMAGIMARY PART OF R_K
       PERT_INI(6)=RPRIME_IM!IMAGINARY PART OF R_K PRIME
     !===================================================================== 
       IF(NFOLD>(NIC+1.0D-10)) THEN

          CALL SETUP(COMM,NIC,PERT_INI,END_OF_INFLATION,TOLERANCE,THRESHOLDS,Method="H",MESSAGE=.FALSE.)
          CALL RANGE_INTEGRATE(COMM,CURV_PERT,NFOLD,NFOLD_GOT,Y_GOT,YDERIV_GOT,FLAG=FLAG)
          IF(FLAG/=1) THEN
            PRINT*,"INTEGRATION NOT SUCCESSFUL"
            PRINT*,'COULD EVALUATE TILL N=',NFOLD_GOT,"INSTEAD N=",NFOLD
            IF(ABS(NFOLD_GOT-NFOLD)>1.0d-10) STOP "EXITING"
          ENDIF
          PERT_INI(:)=Y_GOT(:)
       ENDIF
       
       
       !=============   DEFINTION OF RK3 AND RK3' ==================================== 


       CURV_K2=CMPLX(PERT_INI(3),PERT_INI(5),KIND=WP)

       
       CURV_PRIME_K2=CMPLX(PERT_INI(4),PERT_INI(6),KIND=WP)

       
       CALL COLLECT_GARBAGE(COMM)      

      !=========  At nfold background quantities shall not chage due to change of mode ========
    end if integration_for_scalene
    
    cutoff = DEXP(-TINYVALUE*(MODEK1+MODEK2+MODEK3)/(3.0d0*AI*DEXP(NFOLD)*HUBB))



      IF(TERM==1) THEN 
         INTEGRAND=AI*DEXP(NFOLD)*AI*DEXP(NFOLD)*AI*DEXP(NFOLD)&
	 &*EPSILON1*EPSILON1*HUBB
         FACTOR=2.0d0*((DCONJG(CURV_K1)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K3))&
	&+(DCONJG(CURV_K2)*DCONJG(CURV_PRIME_K3)*DCONJG(CURV_PRIME_K1))&
	&+(DCONJG(CURV_K3)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K2)))*cutoff
        ELSEIF (TERM==2) THEN 
         INTEGRAND=AI*DEXP(NFOLD)*EPSILON1*EPSILON1/HUBB
         FACTOR=-2.0d0*(DCONJG(CURV_K1)*DCONJG(CURV_K2)*DCONJG(CURV_K3))*(k1dotk2+k2dotk3+k3dotk1)*cutoff
	elseif (term==3) then
         INTEGRAND=AI*DEXP(NFOLD)*AI*DEXP(NFOLD)*AI*DEXP(NFOLD)&
	 &*EPSILON1*EPSILON1*HUBB
	 Factor=-2.0d0*(((k1dotk2/modek2**2.0d0)*(DCONJG(CURV_K1)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K3)))&
	 &+((k1dotk2/modek1**2.0d0)*(DCONJG(CURV_K2)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K3)))&
	 &+((k2dotk3/modek3**2.0d0)*(DCONJG(CURV_K2)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K3)))&
	 &+((k2dotk3/modek2**2.0d0)*(DCONJG(CURV_K3)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K2)))&
	 &+((k3dotk1/modek3**2.0d0)*(DCONJG(CURV_K1)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K3)))&
	 &+((k3dotk1/modek1**2.0d0)*(DCONJG(CURV_K3)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K1))))*cutoff
        ELSEIF (TERM==47) THEN 
         INTEGRAND=AI*DEXP(NFOLD)*AI*DEXP(NFOLD)*EPSILON_ETAPRIME    
         FACTOR=(DCONJG(CURV_K1)*DCONJG(CURV_K2)*DCONJG(CURV_PRIME_K3))&
	&+(DCONJG(CURV_K2)*DCONJG(CURV_K3)*DCONJG(CURV_PRIME_K1))&
	&+(DCONJG(CURV_K3)*DCONJG(CURV_K1)*DCONJG(CURV_PRIME_K2))
        ELSEIF (TERM==5) THEN 
         INTEGRAND=AI*DEXP(NFOLD)*AI*DEXP(NFOLD)*AI*DEXP(NFOLD) &
	&*EPSILON1*EPSILON1*EPSILON1*HUBB
	 Factor=0.5d0*(((k1dotk2/modek2**2.0d0)*(DCONJG(CURV_K1)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K3)))&
	 &+((k1dotk2/modek1**2.0d0)*(DCONJG(CURV_K2)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K3)))&
	 &+((k2dotk3/modek3**2.0d0)*(DCONJG(CURV_K2)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K3)))&
	 &+((k2dotk3/modek2**2.0d0)*(DCONJG(CURV_K3)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K2)))&
	 &+((k3dotk1/modek3**2.0d0)*(DCONJG(CURV_K1)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K3)))&
	 &+((k3dotk1/modek1**2.0d0)*(DCONJG(CURV_K3)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K1))))*cutoff
	ELSEIF (TERM==6) THEN
	 INTEGRAND=AI*DEXP(NFOLD)*AI*DEXP(NFOLD)*AI*DEXP(NFOLD) &
	&*EPSILON1*EPSILON1*EPSILON1*HUBB
	 FACTOR=0.5d0*((((modek1*modek1)*k2dotk3/(modek2*modek2*modek3*modek3))&
	 &*(DCONJG(CURV_K1)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K3)))&
	 &+(((modek2*modek2)*k3dotk1/(modek1*modek1*modek3*modek3))&
	 &*(DCONJG(CURV_K2)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K3)))&
	 &+(((modek3*modek3)*k1dotk2/(modek1*modek1*modek2*modek2))&
	 &*(DCONJG(CURV_K3)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K1))))*cutoff
	 ELSE 
	INTEGRAND=0.0D0
	FACTOR=0.0D0
      END IF
      IF(TERM==0) THEN 
	INTEGRAND1=AI*DEXP(NFOLD)*AI*DEXP(NFOLD)*AI*DEXP(NFOLD)&
	 &*EPSILON1*EPSILON1*HUBB*2.0d0&
	&*((DCONJG(CURV_K1)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K3))&
	&+(DCONJG(CURV_K2)*DCONJG(CURV_PRIME_K3)*DCONJG(CURV_PRIME_K1))&
	&+(DCONJG(CURV_K3)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K2)))*cutoff

        INTEGRAND2=AI*DEXP(NFOLD)*EPSILON1*(EPSILON1/HUBB)*(-2.0d0)*&
	&(DCONJG(CURV_K1)*DCONJG(CURV_K2)*DCONJG(CURV_K3))*(k1dotk2+k2dotk3+k3dotk1)*cutoff	

        INTEGRAND3=AI*DEXP(NFOLD)*AI*DEXP(NFOLD)*AI*DEXP(NFOLD)&
	 &*EPSILON1*EPSILON1*HUBB*(-2.0d0)&
	 &*(((k1dotk2/modek2**2.0d0)*(DCONJG(CURV_K1)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K3)))&
	 &+((k1dotk2/modek1**2.0d0)*(DCONJG(CURV_K2)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K3)))&
	 &+((k2dotk3/modek3**2.0d0)*(DCONJG(CURV_K2)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K3)))&
	 &+((k2dotk3/modek2**2.0d0)*(DCONJG(CURV_K3)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K2)))&
	 &+((k3dotk1/modek3**2.0d0)*(DCONJG(CURV_K1)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K3)))&
	 &+((k3dotk1/modek1**2.0d0)*(DCONJG(CURV_K3)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K1))))*cutoff

        INTEGRAND4=AI*DEXP(NFOLD)*AI*DEXP(NFOLD)*EPSILON_ETAPRIME&
	 &*((DCONJG(CURV_K1)*DCONJG(CURV_K2)*DCONJG(CURV_PRIME_K3))&
	 &+(DCONJG(CURV_K2)*DCONJG(CURV_K3)*DCONJG(CURV_PRIME_K1))&
	 &+(DCONJG(CURV_K3)*DCONJG(CURV_K1)*DCONJG(CURV_PRIME_K2)))           
  
	INTEGRAND5=AI*DEXP(NFOLD)*AI*DEXP(NFOLD)*AI*DEXP(NFOLD) &
	 &*EPSILON1*EPSILON1*EPSILON1*HUBB*0.5d0*(((k1dotk2/modek2**2.0d0)*&
	 &(DCONJG(CURV_K1)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K3)))&
	 &+((k1dotk2/modek1**2.0d0)*(DCONJG(CURV_K2)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K3)))&
	 &+((k2dotk3/modek3**2.0d0)*(DCONJG(CURV_K2)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K3)))&
	 &+((k2dotk3/modek2**2.0d0)*(DCONJG(CURV_K3)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K2)))&
	 &+((k3dotk1/modek3**2.0d0)*(DCONJG(CURV_K1)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K3)))&
	 &+((k3dotk1/modek1**2.0d0)*(DCONJG(CURV_K3)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K1))))*cutoff
        
        INTEGRAND6=AI*DEXP(NFOLD)*AI*DEXP(NFOLD)*AI*DEXP(NFOLD)&
       	 &*EPSILON1*EPSILON1*EPSILON1*HUBB*0.5d0*((((modek1*modek1)*k2dotk3/(modek2*modek2*modek3*modek3))&
	 &*(DCONJG(CURV_K1)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K3)))&
	 &+(((modek2*modek2)*k3dotk1/(modek1*modek1*modek3*modek3))&
	 &*(DCONJG(CURV_K2)*DCONJG(CURV_PRIME_K1)*DCONJG(CURV_PRIME_K3)))&
	 &+(((modek3*modek3)*k1dotk2/(modek1*modek1*modek2*modek2))&
	 &*(DCONJG(CURV_K3)*DCONJG(CURV_PRIME_K2)*DCONJG(CURV_PRIME_K1))))*cutoff
      end if   
!=======================  changed =======================================
       
       
       IF(RL_INTEGRATION)THEN
          
          INTEGRAND=INTEGRAND*REAL(FACTOR)
          if(TERM==0)then
          INTEGRAND=REAL(INTEGRAND1+INTEGRAND2+INTEGRAND3+INTEGRAND4+INTEGRAND5+INTEGRAND6)
          end if
       END IF
           
       IF(CMPLX_INTEGRATION)THEN
          
          INTEGRAND=INTEGRAND*AIMAG(FACTOR)
          if(TERM==0)then
          INTEGRAND=AIMAG(INTEGRAND1+INTEGRAND2+INTEGRAND3+INTEGRAND4+INTEGRAND5+INTEGRAND6)
          end if
          
       END IF
       IF(RL_INTEGRATION.AND.CMPLX_INTEGRATION) STOP 'CHECK THE REAL AND COMPLEX LOGICAL SWITCH'
       
       

       IF(SAVE_RK_NEND)THEN
          RK1_NEND=CURV_K1
          RK1_PRIME_NEND=CURV_PRIME_K1
	  if(equilateral) then
	  RK2_NEND=RK1_NEND
	  RK2_PRIME_NEND=RK1_PRIME_NEND
	  RK3_NEND=RK1_NEND
	  RK3_PRIME_NEND=RK1_PRIME_NEND
	  end if 
	  if(isosceles) then
	  RK2_NEND=RK1_NEND
	  RK2_PRIME_NEND=RK1_PRIME_NEND
	  RK3_NEND=CURV_K3
	  RK3_PRIME_NEND=CURV_PRIME_K3
	  end if 
          if(scalene)then
	  RK2_NEND=CURV_K2
	  RK2_PRIME_NEND=CURV_PRIME_K2
	  RK3_NEND=CURV_K3
	  RK3_PRIME_NEND=CURV_PRIME_K3
          endif
       END IF

 
     END FUNCTION INTEGRAND
