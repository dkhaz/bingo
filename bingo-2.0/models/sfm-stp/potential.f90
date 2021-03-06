!DEFINE YOUR POTENTIAL HERE
     MODULE POT2DER
       USE RKSUITE_90_PREC

       REAL(WP)::V_PHI_PHI
     END MODULE POT2DER


     SUBROUTINE POTENTIAL_POTENTIALPRIME(POTENTIAL_USED,POTENTIALPRIME,PHI_X)
       USE RKSUITE_90_PREC
       USE THEORY_PARAM
       USE PARAMETERS
       USE POT2DER
       IMPLICIT NONE

       REAL(WP),INTENT(IN)::PHI_X
       REAL(WP),INTENT(OUT)::POTENTIAL_USED,POTENTIALPRIME
       REAL(WP)::V00,ALPHA,PHI0,DELTA_PHI 
       POTENTIAL_USED=0.0_WP
       POTENTIALPRIME=0.0_WP
       V_PHI_PHI=0.0_WP

       V00=PARAM_1
       ALPHA=PARAM_2
       PHI0=PARAM_3
       DELTA_PHI=PARAM_4
       
       POTENTIAL_USED=V00*(1.0_WP-((PHI_X/MU)**PP))*(1.0_WP+ALPHA*TANH((PHI_X-PHI0)/DELTA_PHI))
       POTENTIALPRIME=-(V00/MU)*PP*((PHI_X/MU)**(PP-1.0_WP))*(1.0_WP+ALPHA*TANH((PHI_X-PHI0)/DELTA_PHI))
       POTENTIALPRIME=POTENTIALPRIME+(V00*(1.0_WP-((PHI_X/MU)**PP))*(ALPHA/DELTA_PHI)/(COSH((PHI_X-PHI0)/DELTA_PHI)&
       &*COSH((PHI_X-PHI0)/DELTA_PHI)))
  
       V_PHI_PHI=-(V00*PP*(PP-1.0D0)/(MU*MU))*(1.0_WP+ALPHA*TANH((PHI_X-PHI0)/DELTA_PHI))*((PHI_X/MU)**(PP-2.0D0))
  
       V_PHI_PHI=V_PHI_PHI-2.0D0*(V00*PP*ALPHA/(MU*DELTA_PHI))*((PHI_X/MU)**(PP-1.0D0))&
       &*(1.0D0/(COSH((PHI_X-PHI0)/DELTA_PHI)*COSH((PHI_X-PHI0)/DELTA_PHI)))
  
       V_PHI_PHI=V_PHI_PHI-2.0D0*(ALPHA*V00/(DELTA_PHI*DELTA_PHI))&
       &*(1.0D0-(PHI_X/MU)**PP)*(TANH((PHI_X-PHI0)/DELTA_PHI)/(COSH((PHI_X-PHI0)/DELTA_PHI)*COSH((PHI_X-PHI0)/DELTA_PHI)))
    
     
     END SUBROUTINE POTENTIAL_POTENTIALPRIME
