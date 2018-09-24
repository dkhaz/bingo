#ifort compiler (11 or later)
#Assuming mpif90 points to ifort
   F90C     = mpif90
   FFLAGS   = -O3 -openmp -fpp -DMPI
#ifort compiler (11 or later)

#   F90C     = ifort
#   FFLAGS   = -O3 -openmp -fpp 
#  -warn noerrors
#============== Others (auto) ================
#   F90C     = mpif90
#   FFLAGS   = -O3 $(shell mpif90 --showme:compile) $(shell mpif90 --showme:link)
#=============================================
#   F90C     = gfortran
#   FFLAGS   = -O3 -cpp
#=============================================

#==============================================
#------------- Models--------------------------
 Model = qp-stp
# Model =  plaw
# Model = sfm-stp
# Model = amm
# Model = pi
# Model = Sm
#==============================================


Modelfolder  = models/$(Model)

Theory = theory_param.o

Parameters = parameters.o

Potential = potential.o

Modelfiles =  $(Theory) $(Parameters) $(Potential)

NDSOLVER = rksuite.o

INTIGRATOR = d1mach.o dqpsrt.o dqk21.o dqelg.o dqagse.o

OBJFILES = inifile.o store.o interfacing.o eqn_bgrnd.o bgrnd.o \
	eqn_pert.o find_nic.o find_nshs.o initial_cond.o dotprod.o integrand_N.o $(INTIGRATOR) driver.o

OBJFILES3D = inifile.o store.o interfacing.o eqn_bgrnd.o bgrnd.o \
	eqn_pert.o find_nic.o find_nshs.o initial_cond.o dotprod.o integrand_N.o $(INTIGRATOR) 3driver.o



export FFLAGS
export F90C


%.o: %.f90
	$(F90C) $(FFLAGS) -c $*.f90
%.o: %.f
	$(F90C) $(FFLAGS) -c $*.f
all :	bingo2

gather :
	$(F90C) gather.f90 -o gatherall.out

bingo2:	$(NDSOLVER) $(Modelfolder)/$(Theory) $(Modelfolder)/$(Parameters) $(Modelfolder)/$(Potential) $(OBJFILES)
	$(F90C) $(FFLAGS) $(NDSOLVER) $(Modelfiles) $(OBJFILES) -o bingo2.out 

bingo3d:$(NDSOLVER) $(Modelfolder)/$(Theory) $(Modelfolder)/$(Parameters) $(Modelfolder)/$(Potential) $(OBJFILES)
	$(F90C) $(FFLAGS) $(NDSOLVER) $(Modelfiles) $(OBJFILES3D) -o bingo2.out 


run :
	time ./bingo2.out $(Modelfolder)/fnlparam.ini 
prun:
	mpirun -np 2 time ./bingo2.out $(Modelfolder)/fnlparam.ini
figs:
	gnuplot > load 'plots.p'
view :
	gv allplots.ps


clean:
	rm -f *.o *.mod *.d *.pc *.obj fort.* *.out *.dat *.p~ *.f90~ Makefile~ *.ini~ plots/* ../core

