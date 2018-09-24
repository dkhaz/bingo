! The main driver. 

PROGRAM EVALUATION
  
  USE INIFILE 
  USE RKSUITE_90_PREC
  USE THEORY_PARAM
  USE PARAMETERS
  USE INTERFACING
  USE BACKGROUND_PATH
  USE STOREK
  USE FNLTERM
  USE MODESTORE
  USE RKSUITE_90_PREC
  IMPLICIT NONE
#ifdef MPI
  include 'mpif.h'
#endif  
  REAL(WP)::INTATEND,K1EXP,K2EXP,KEXP,INCRMNT,abserr,G7,PERM7,REALGTOT,PERMUTATION
  REAL(WP)::INTSTART,INTEND,epsrel,epsabs,ANS_RL,ANS_CMPLX,ANS,onetwo,twothree,threeone,onebythree,twobythree,harvest
  REAL(WP),EXTERNAL::INTEGRAND
  Real(wp)::ini_pot,ini_pot_prime,ini_hub,smallestmode,largestmode
  real(wp),allocatable,dimension(:)::alist,blist,rlist,elist,modes,arrayktwobykone
  REAL(WP)::POWERSPECTRA1,POWERSPECTRA2,POWERSPECTRA3
  integer,allocatable,dimension(:)::iord 
  integer,parameter::limit_subintervals=20000 
  INTEGER::IFAIL,kpoints,neval,last,i,j,modeelements,jend,jstart,ISTART,IEND,INTERVAL,kk
  character(LEN=Ini_max_string_len)::InputFile,GFILENAME,FILENUM
  COMPLEX(WP)::ID,F_NL,IMAGINARY,G123
  logical:: bad,initial_slow_roll,calc_fnl,squeezed,triangle_form,dummy_print,compute_triangle
#ifdef MPI
  INTEGER ( KIND = 4 ) ::ERROR_FLAG,P_NUM,MY_ID,PFILENUM
  character(LEN=Ini_max_string_len)::FNLFILENAME,numstr
  logical::seed_first_time=.true.
#endif
  !=================== Input file   =========================
  InputFile = ''
  
  if (iargc() /= 0)  call getarg(1,InputFile)
  if (InputFile == '') stop 'No parameter input file'

  call Ini_Open(InputFile, 1, bad, .false.)
  if (bad) stop 'Error opening parameter file'
  Ini_fail_on_not_found = .false.
   
  !=========================================================
  IMAGINARY=(0.0_WP,1.0_WP)
  F_NL=(0.0D0,0.0D0)    

  
  PRINT*,"===================================================================="
  PRINT*,'THIS PROGRAMME EVALUATES THE POWER SPECTRUM AND THE BISPECTRUM &
       &  FOR A GIVEN FIELD AND POTENTIAL'
  PRINT*,"===================================================================="
  PRINT*,'THE FIELD USED HERE : CANONICAL'
  PRINT*,"===================================================================="
  !========  Reading from file ===============================
  param_1=Ini_Read_double('param1',0.0d0)
  param_2=Ini_Read_double('param2',0.0d0)
  param_3=Ini_Read_double('param3',0.0d0)
  param_4=Ini_Read_double('param4',0.0d0)
  
  initialphi=Ini_Read_double('phi_i') 
  initial_slow_roll=Ini_Read_Logical('start_with_slow_roll')
  
  if (initial_slow_roll) then
     !============== Taking it to the attractor=====================================
     
     call potential_potentialprime(ini_pot,ini_pot_prime,initialphi)
     
     !=============  slow roll determines the phi_dot ====================
     
     ini_hub=dsqrt(ini_pot/3.0_wp)
     INITIALPHIDOT=-ini_pot_prime/(3.0_wp*ini_hub)
     
     else
        INITIALPHIDOT=Ini_Read_double('phi_dot_i')
        
  end if
  

  pivotk=Ini_Read_double('pivot_scale',0.05d0)
  N_pivot=Ini_Read_double('Npivot',50.0d0)      
  
  force_ai=Ini_Read_Logical('force_aini')
  if(force_ai) AFI=Ini_Read_double('ainitial')
  multiphase=Ini_Read_Logical('multi_phase')
  if(multiphase) mpe=Ini_Read_double('expected_multi_phase_end')
  acc=Ini_Read_double('accuracy',0.0d0)
  !===========================================================
  CALL BACKGROUND_VALUE
  PRINT*,'THE BACKGROUND OF THE MODEL HAS BEEN EVALUATED'
  
  !===========================================================================
  
  epsrel=(1.0D-5)*(10.0**(-acc))
  epsabs=0.0d0
  equilateral=Ini_Read_Logical('Equilateral')
  isosceles=Ini_Read_Logical('Isosceles')
  scalene=Ini_Read_Logical('Scalene')
  squeezed=Ini_Read_Logical('Squeezed')

  compute_triangle=.false.
  
  if(scalene) compute_triangle=Ini_Read_Logical('Compute_traiangle',.False.)
  
  if(squeezed) isosceles=.True.
  
#ifdef MPI
  if(.not.(scalene)) stop'MPI can be used only in scalene case'
#endif
#ifdef MPI
  call MPI_Init ( error_flag )
   if (error_flag/=MPI_SUCCESS) stop 'MPI fail: rank'
  call MPI_Comm_size ( MPI_COMM_WORLD, p_num, error_flag )
  call MPI_Comm_rank ( MPI_COMM_WORLD, my_id, error_flag )
  PFILENUM=MY_ID+1
  write(numstr,'(i10)') PFILENUM
  
  FNLFILENAME=trim("F_nl")//"_"//trim(adjustl(numstr))//".txt"
#endif
  ! ======= INTEGRATION STARTS HERE =========================================
  
  K1EXP=Ini_Read_Double('logki',-5.0d0)
  K2EXP=Ini_Read_Double('logkf',-1.0d0)
  kpoints=Ini_Read_Int('num_k',500)
  INCRMNT=(K2EXP-K1EXP)/real(kpoints)
  calc_fnl=Ini_Read_Logical('calcfnl')
  term=Ini_Read_Int('Term',47)
  !============ changed =============================  
  

  if((equilateral.and.isosceles).or.(equilateral.and.scalene)&
	&.or.(isosceles.and.scalene)) stop "Triangle Not properly defined : fnlparams.ini"
  
  if(.not.(equilateral)) then
     modek3=Ini_Read_Double('Fixedmodexpo',-5.0d0)
     modek3=(10.0d0**modek3)
  end if
  if((.not.(equilateral)).and.(.not.(calc_fnl))) stop "To calculate only &
       & power spectrum set Equilateral= T"
  allocate(modes(kpoints))
  do i=1,kpoints
     modes(i)=k1exp+(k2exp-k1exp)*dble(i-1)/dble(kpoints-1)
  enddo
  
  modes=(10.0d0)**modes
  
  IF(calc_fnl.and.(TERM/=1).and.(TERM/=2).and.(TERM/=3)&
       &.and.(TERM/=47).and.(TERM/=5).and.(TERM/=6).and.(TERM/=0)) stop "Specify the term"
  
  NIC_COND=Ini_Read_double('Nicond',100.0d0)

  !========== For calculating the largest mode that is sub-Hubble at the beginning of the evolution =================
  if(squeezed) MODEK3=NIC_COND*ai*dexp(ef(2))*Hubble(2)
  !========== For calculating the largest mode that is sub-Hubble at the beginning of the evolution =================
  
    
  !============ changed =============================  
  if(equilateral) then 
     write(FILENUM,'(i10)') TERM
     GFILENAME=trim("k6-G_")//trim(adjustl(FILENUM))//".txt"
     OPEN(UNIT=901,FILE='plots/PS.txt',FORM='FORMATTED',STATUS='REPLACE')
     OPEN(UNIT=902,FILE='plots/'//GFILENAME,FORM='FORMATTED',STATUS='REPLACE')
  end if
  
#ifndef MPI  
  OPEN(UNIT=903,FILE='plots/F_NL-RE.txt',FORM='FORMATTED',STATUS='REPLACE')
  ISTART=1
  IEND=KPOINTS
  INTERVAL=1
#endif
#ifdef MPI
  OPEN(UNIT=903+MY_ID,FILE='plots/'//FNLFILENAME,FORM='FORMATTED',STATUS='REPLACE')
  ISTART=MY_ID+1
  IEND=KPOINTS-P_NUM+MY_ID+1
  INTERVAL=P_NUM
#endif 


  DO I=ISTART,IEND,INTERVAL
     
     MODEK2=MODES(I)
     
     DO J=1,kpoints
     
        MODEK1=MODES(J)        

     Do kk=1,kpoints

         modek3=MODES(kk) 
       
        acc_check:if(largestmode>2.0d-2) then 
           epsrel=(1.0D-5)*(10.0**(-acc))
           acc_check_2:if(acc>2) then
              if(NIC_COND<250.0d0) NIC_COND=250.0D0*(dble(acc)/2.0d0)
           else
              if(NIC_COND<250.0d0) NIC_COND=250.0D0
           endif acc_check_2
          !========== For calculating the largest mode that is sub-Hubble at the beginning of the evolution =================
	  if(squeezed) MODEK3=NIC_COND*ai*dexp(ef(2))*Hubble(2)
	  !========== For calculating the largest mode that is sub-Hubble at the beginning of the evolution =================
  
        endif acc_check

        if(equilateral.or.isosceles) modek2=modek1
	if(squeezed.and.(modek1<=modek3)) stop "Squeezed Condition not satisfied"   
	if(equilateral) modek3=modek1
        nic_nshs_mode:if(.not.(equilateral)) then    
           smallestmode=min(min(modek1,modek2),modek3) 
           largestmode=max(max(modek1,modek2),modek3)
        else
           smallestmode=modek1 
           largestmode=modek1
        end if nic_nshs_mode
       
        
        call dotproduct(onetwo,twothree,threeone,modek1,modek2,modek3)
        k1dotk2=onetwo
        k2dotk3=twothree
        k3dotk1=threeone
        
        CALL FIND_NIC(largestmode)
        
        CALL FIND_SHS(largestmode)
        
        
        INTSTART=NIC
        INTEND=NSHS
	triangle_form=.False.
        if(modek1+modek2>modek3.and.modek2+modek3>modek1.and.modek3+modek1>modek2) triangle_form=.true.
        calc_fnl_check1:if(calc_fnl.and.triangle_form) then
#ifndef MPI
	   PRINT*,"CARRYING OUT INTEGRATION FOR THE MODE =",MODEK1,MODEK2,MODEK3
#endif
           F_NL=(0.0D0,0.0D0)
           
           SAVE_RK_NEND=.FALSE.
           RL_INTEGRATION=.TRUE.
           CMPLX_INTEGRATION=.FALSE.
           
           ANS_RL=0.0D0
           ANS_CMPLX=0.0D0
           
           allocate(alist(limit_subintervals),blist(limit_subintervals),rlist(limit_subintervals) &
                &,elist(limit_subintervals),iord(limit_subintervals)) 
           call dqagse(INTEGRAND,INTSTART,INTEND,epsabs,epsrel,limit_subintervals,&
                &ans,abserr,neval,ifail,alist,blist,rlist,elist,iord,last)
           
           deallocate(alist,blist,rlist,elist,iord)
           ANS_RL=ans
     
           
           
           RL_INTEGRATION=.FALSE.
           CMPLX_INTEGRATION=.TRUE.

           allocate(alist(limit_subintervals),blist(limit_subintervals),rlist(limit_subintervals) &
                &,elist(limit_subintervals),iord(limit_subintervals)) 
           
           call dqagse(INTEGRAND,INTSTART,INTEND,epsabs,epsrel,limit_subintervals,&
                &ans,abserr,neval,ifail,alist,blist,rlist,elist,iord,last)
           
           deallocate(alist,blist,rlist,elist,iord)
           
           ANS_CMPLX=ans
           
           ID=IMAGINARY*CMPLX(ANS_RL,ANS_CMPLX,KIND=WP)
           
        endif calc_fnl_check1

        IF(triangle_form) THEN
        SAVE_RK_NEND=.TRUE.
        INTATEND=INTEGRAND(INTEND)
        
        POWERSPECTRA1=2.0_WP*(ABS(RK1_NEND)*ABS(RK1_NEND))*((MODEK1**3.0_WP)/(TWOPI_D**2.0_WP))
        POWERSPECTRA2=2.0_WP*(ABS(RK2_NEND)*ABS(RK2_NEND))*((MODEK2**3.0_WP)/(TWOPI_D**2.0_WP))
        POWERSPECTRA3=2.0_WP*(ABS(RK3_NEND)*ABS(RK3_NEND))*((MODEK3**3.0_WP)/(TWOPI_D**2.0_WP))
        
        PERMUTATION=((MODEK1*MODEK1*MODEK1)*POWERSPECTRA2*POWERSPECTRA3)&
             &+((MODEK2*MODEK2*MODEK2)*POWERSPECTRA3*POWERSPECTRA1)&
             &+((MODEK3*MODEK3*MODEK3)*POWERSPECTRA1*POWERSPECTRA2)
        end if

        calc_fnl_check:IF(calc_fnl.and.triangle_form) THEN 
           G123=((RK1_NEND*RK2_NEND*RK3_NEND)*ID)+DCONJG((RK1_NEND*RK2_NEND*RK3_NEND)*ID)
           
           IF(TERM==47.or.TERM==0) THEN 
              PERM7=ABS(RK2_NEND)*ABS(RK2_NEND)*ABS(RK3_NEND)*ABS(RK3_NEND)&
                   &+ABS(RK3_NEND)*ABS(RK3_NEND)*ABS(RK1_NEND)*ABS(RK1_NEND)&
                   &+ABS(RK1_NEND)*ABS(RK1_NEND)*ABS(RK2_NEND)*ABS(RK2_NEND)
              
              G7=PERM7*EPS2*0.5d0
           ELSE
              G7=0.0D0
           END IF
           
           REALGTOT=(dble(real(G123))+G7)
           if(equilateral) WRITE(902,*)MODEK1,REALGTOT*MODEK1**6.0D0  
           
           
           F_NL=(-10.0D0/3.0D0)*((TWOPI_D)**(-4.0d0))*(REALGTOT*MODEK1*MODEK2 &
                &*MODEK3*MODEK1*MODEK2*MODEK3*MODEK1*MODEK2*MODEK3)
           F_NL=F_NL/PERMUTATION
#ifndef MPI 
           PRINT*,'F_NL=',DBLE(F_NL)
#endif
           if(equilateral.or.isosceles) WRITE(903,*)MODEK1,DBLE(F_NL)
#ifndef MPI
          if(scalene) WRITE(903,'(999E15.5)')modek1,modek2,modek3,DBLE(F_NL)
#endif
#ifdef MPI
          if(scalene) WRITE(903+MY_ID,'(999E15.5)')modek1,modek2,modek3,DBLE(F_NL)
#endif

        END IF calc_fnl_check
     IF(calc_fnl.and.(.not.(triangle_form))) THEN 
#ifndef MPI
          if(scalene) WRITE(903,'(999E15.5)')modek1,modek2,modek3,0
#endif
#ifdef MPI
          if(scalene) WRITE(903+MY_ID,'(999E15.5)')modek1,modek2,modek3,0
#endif

      END IF		
        if(equilateral) WRITE(901,*)MODEK1,POWERSPECTRA1
        RL_INTEGRATION=.FALSE.
        CMPLX_INTEGRATION=.FALSE.
     END DO
     if(.not.(Compute_triangle)) Exit 

  END DO

  END DO
#ifdef MPI
  close(903+MY_ID)
  call MPI_Finalize (error_flag )
  
#endif
  if(equilateral) CLOSE(901)
  if(equilateral) CLOSE(902)
#ifndef MPI
  CLOSE(903)
#endif
  !!==========================================================================
  PRINT*,'THE PROGRAME IS COMPLETED'
  
END PROGRAM EVALUATION
