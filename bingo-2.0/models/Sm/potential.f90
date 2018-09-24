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
       REAL(WP)::V00,plus,minus,tanhfucn,sechfunc,yp,fs,aplus,aminus,phi0
 
       POTENTIAL_USED=0.0_WP
       POTENTIALPRIME=0.0_WP
       V_PHI_PHI=0.0_WP

       
       V00=PARAM_1
       aplus=PARAM_2
       aminus=PARAM_3
       phi0=PARAM_4
       
       fs=1d4/(INITIALPHI-phi0)
       
       yp=phi_x-phi0
       plus=aplus+aminus
       minus=aplus-aminus

       
       sechfunc=1.0d0/(dcosh(fs*yp))
       tanhfucn=dtanh(fs*yp)
 
       potential_used=V00+(0.5d0*plus*yp)+(0.5d0*yp*tanhfucn*minus)

       
       potentialprime=(0.5d0*plus)+(0.5d0*tanhfucn*minus)+0.5d0*fs*yp*sechfunc*sechfunc*minus
       

       v_phi_phi=(fs*minus*sechfunc*sechfunc)-fs*fs*yp*sechfunc*sechfunc*tanhfucn*minus
           
     
     END SUBROUTINE POTENTIAL_POTENTIALPRIME
