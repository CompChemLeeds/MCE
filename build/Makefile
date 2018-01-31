#Makefile for use on chmlinxx computers.
#Made with the aid of the fmkmf utility (Copyright (c) The University of Edinburgh) 

F90=gfortran
FFLAGS=-fopenmp -O3
FLINKS=-lm -Warray-bounds
VPATH=../src

# Link step

all: Main avrgpops timehist subavrg interpolate integrator

Main:BLAS.o zgesv_with_dpndnss.o zheev_with_dpndnss.o randgen.o dcerf.o vars.o alarrays.o sb.o vp.o hp.o fp.o mp.o iv.o cp.o hh.o outputs.o redirect.o Ham.o Chks.o bsetalter.o readpars.o derivsMCE.o propMCE.o bsetgen.o MainMCE.o 
	 $(F90) $(FFLAGS) -o MCE.exe BLAS.o zgesv_with_dpndnss.o zheev_with_dpndnss.o randgen.o dcerf.o vars.o alarrays.o sb.o vp.o hp.o fp.o mp.o iv.o cp.o hh.o outputs.o redirect.o Ham.o Chks.o bsetalter.o readpars.o derivsMCE.o propMCE.o bsetgen.o MainMCE.o $(FLINKS)

avrgpops:avrgpops.f90
	 $(F90) -O3 -o avrgpops.exe $<

subavrg:subavrg.f90
	 $(F90) -O3 -o subavrg.exe $<

timehist:timehist.f90
	 $(F90) -O3 -o timehist.exe $<

interpolate:neville.o interpolate.o
	 $(F90) -O3 -o interpolate.exe neville.o interpolate.o

integrator:integrator.f90
	 $(F90) -O3 -o integrator.exe $<

# Compile steps

BLAS.o:BLAS.f
	 $(F90) $(FFLAGS) -c $<

zgesv_with_dpndnss.o:zgesv_with_dpndnss.f
	 $(F90) $(FFLAGS) -c $<

zheev_with_dpndnss.o:zheev_with_dpndnss.f
	 $(F90) $(FFLAGS) -c $<

randgen.o:randgen.f
	 $(F90) $(FFLAGS) -c $<

dcerf.o:dcerf.f90  
	 $(F90) $(FFLAGS) -c $<

vars.o:vars.f90  
	 $(F90) $(FFLAGS) -c $<

alarrays.o:alarrays.f90 vars.o 
	 $(F90) $(FFLAGS) -c $<

sb.o:sb.f90 vars.o 
	 $(F90) $(FFLAGS) -c $<

vp.o:vp.f90 vars.o 
	 $(F90) $(FFLAGS) -c $<

hp.o:hp.f90 vars.o 
	 $(F90) $(FFLAGS) -c $<

fp.o:fp.f90 vars.o 
	 $(F90) $(FFLAGS) -c $<

mp.o:mp.f90 vars.o 
	 $(F90) $(FFLAGS) -c $<

iv.o:iv.f90 vars.o 
	 $(F90) $(FFLAGS) -c $<

cp.o:cp.f90 dcerf.o vars.o 
	 $(F90) $(FFLAGS) -c $<

hh.o:hh.f90 vars.o 
	 $(F90) $(FFLAGS) -c $<

outputs.o:outputs.f90 vars.o 
	 $(F90) $(FFLAGS) -c $<

redirect.o:redirect.f90 vars.o sb.o vp.o hp.o fp.o mp.o iv.o cp.o hh.o 
	 $(F90) $(FFLAGS) -c $<

Ham.o:Ham.f90 vars.o alarrays.o redirect.o 
	 $(F90) $(FFLAGS) -c $<

Chks.o:Chks.f90 vars.o Ham.o redirect.o
	 $(F90) $(FFLAGS) -c $<

bsetalter.o:bsetalter.f90 vars.o Ham.o alarrays.o outputs.o 
	 $(F90) $(FFLAGS) -c $<

readpars.o:readpars.f90 vars.o alarrays.o Ham.o redirect.o 
	 $(F90) $(FFLAGS) -c $<

derivsMCE.o:derivsMCE.f90 vars.o Ham.o alarrays.o outputs.o redirect.o readpars.o 
	 $(F90) $(FFLAGS) -c $<

propMCE.o:propMCE.f90 alarrays.o vars.o derivsMCE.o
	 $(F90) $(FFLAGS) -c $<

bsetgen.o:bsetgen.f90 vars.o Ham.o alarrays.o outputs.o Chks.o redirect.o propMCE.o 
	 $(F90) $(FFLAGS) -c $<

MainMCE.o:MainMCE.f90 bsetgen.o bsetalter.o Ham.o vars.o readpars.o outputs.o alarrays.o Chks.o propMCE.o redirect.o 
	 $(F90) $(FFLAGS) -c $<

neville.o:neville.f90
	 $(F90) -O2 -c $<

interpolate.o:interpolate.f90 neville.o
	 $(F90) -O2 -c $<

clean:
	rm -rf f_{files,modd}* *.o *.mod *.exe *~ ../run/*.exe ../run/*~ ../src/*~

