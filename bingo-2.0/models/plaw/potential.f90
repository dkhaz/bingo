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
       REAL(WP)::V00,qq

       POTENTIAL_USED=0.0_WP
       POTENTIALPRIME=0.0_WP
       V_PHI_PHI=0.0_WP

       
       V00=PARAM_1
       qq=PARAM_2
       
       
       potential_used=v00*dexp(-(dsqrt(2.0_wp/qq))*(phi_x-initialphi))
       potentialprime=-(dsqrt(2.0_wp/qq))*potential_used

       v_phi_phi=(2.0_wp/qq)*potential_used
           
     
     END SUBROUTINE POTENTIAL_POTENTIALPRIME
