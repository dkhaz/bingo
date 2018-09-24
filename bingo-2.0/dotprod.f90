subroutine dotproduct(k1k2,k2k3,k3k1,k1,k2,k3)
 USE RKSUITE_90_PREC

 IMPLICIT NONE

 REAL(WP), INTENT(IN)::k1,k2,k3
 REAL(WP), INTENT(OUT)::k1k2,k2k3,k3k1
 
 k1k2=0.5d0*((k3*k3)-(k1*k1)-(k2*k2))
 k2k3=0.5d0*((k1*k1)-(k2*k2)-(k3*k3))
 k3k1=0.5d0*((k2*k2)-(k3*k3)-(k1*k1))
 

end subroutine dotproduct