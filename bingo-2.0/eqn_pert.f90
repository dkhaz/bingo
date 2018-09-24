! THE CURVATURE PERTURBATION EQUATION
       FUNCTION CURV_PERT(XC,YC)
         USE RKSUITE_90_PREC
         USE PARAMETERS
         USE BACKGROUND_PATH
         USE THEORY_PARAM
         USE ODE_K
         USE INTERFACING, ONLY:POTENTIAL_POTENTIALPRIME
         
         IMPLICIT NONE
         REAL(WP), INTENT(IN) :: XC
         REAL(WP), DIMENSION(:), INTENT(IN) :: YC
         REAL(KIND=WP), DIMENSION(SIZE(YC)) :: CURV_PERT
         REAL(WP)::VC,V_PHIC,HUBC,ZPRIMEC,ZC,AC


         !------------------------BACKGROUND EQUATION-------------------------------
         
         CALL  POTENTIAL_POTENTIALPRIME(VC,V_PHIC,YC(1))

         CURV_PERT(1)=YC(2)
         CURV_PERT(2)=-(3.0_WP-0.5_WP*(YC(2)**2.0_WP))*YC(2)-(6.0_WP-(YC(2)**2.0_WP))*V_PHIC/(2.0_WP*VC)
         !-------------------------------------------------------------------------


         AC=AI*DEXP(XC)


         
         HUBC=VC/(3.0_WP-0.5_WP*(YC(2)**2.0_WP))
         
         HUBC=DSQRT(HUBC)

         ZC=AC*YC(2)
         ZPRIMEC=AC*(YC(2)+CURV_PERT(2))


         !--------------------------CURVATURE PERTURBATION EQUATION----------------------
         CURV_PERT(3)=YC(4)!THIS IS R_K' REAL
         CURV_PERT(4)=-(1.0_WP-(YC(2)**2.0_WP)*0.5_WP+2.0_WP*(ZPRIMEC/ZC))*YC(4)-((MODEK/(AC*HUBC))**2.0_WP)*YC(3)
         CURV_PERT(5)=YC(6)!THIS IS R_K' COMPLEX
         CURV_PERT(6)=-(1.0_WP-(YC(2)**2.0_WP)*0.5_WP+2.0_WP*(ZPRIMEC/ZC))*YC(6)-((MODEK/(AC*HUBC))**2.0_WP)*YC(5)

         !---------------------------------------------------------------------------------


       END FUNCTION CURV_PERT
