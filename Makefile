SFVT_DIR = sfvt
SGM_DIR = sgm
PAMM_DIR = pamm
RTM_DIR = rtm
LSM_DIR = lsm

.PHONY: sfvt_code sgm_code pamm_code rtm_code lsm_code cleanall
CLEANDIRS = $(SFVT_DIR) $(SGM_DIR) $(PAMM_DIR) $(PAMM_DIR)/sfvt $(PAMM_DIR)/bam \
    $(RTM_DIR) $(RTM_DIR)/pts $(LSM_DIR) ./


DEBUG = -fbounds-check -g 
OPT    =-O3

# these three lines should be edited for your system. On systems 
# that do not have separate fortran and c libraries, set NETCDF_FOR and NETCDF_C
# to the same, and set NETCDF_LIB to -lnetcdf (i.e. without the extra f)
#NETCDF_FOR=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.4-mac/
#NETCDF_C=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.1.1-mac/
NETCDF_LIB=-lnetcdff 

NETCDFLIB=-L ${NETCDF_FOR}/lib/  \
          -L ${NETCDF_C}/lib/
NETCDFMOD= ${NETCDF_FOR}/include/


FOR = mpif90 -c  
FOR2 = mpif90  

AR = ar 
RANLIB = ranlib 
OBJ = o

FFLAGS = $(OPT)  $(DEBUG) -w -o 
FFLAGSOMP = -fopenmp-simd $(FFLAGS)
FFLAGS2 =  $(DEBUG) -w -O3 -o 


main.exe	:  main.$(OBJ) variables.$(OBJ) nrtype.$(OBJ) mpi_module.$(OBJ) \
            diagnostics.$(OBJ) \
			 initialisation.$(OBJ) driver_code.$(OBJ) \
			  dynamics.$(OBJ) model_lib.a sfvt_code pamm_code rtm_code lsm_code
	$(FOR2) $(FFLAGSOMP)main.exe main.$(OBJ) variables.$(OBJ) mpi_module.$(OBJ) \
	    diagnostics.$(OBJ) \
		 initialisation.$(OBJ) driver_code.$(OBJ) \
		  dynamics.$(OBJ) \
		  $(SFVT_DIR)/model_lib.a \
		  $(SGM_DIR)/sg_model_lib.a \
		  $(PAMM_DIR)/pmicro_lib.a \
		  $(RTM_DIR)/rtm_lib.a \
		  $(LSM_DIR)/lsm_lib.a \
		  -lm model_lib.a \
		 ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
model_lib.a	:   nrtype.$(OBJ) nr.$(OBJ) nrutil.$(OBJ) locate.$(OBJ) polint.$(OBJ) \
				rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ) zbrent.$(OBJ) \
				hygfx.$(OBJ)  random.$(OBJ) 
	$(AR) rc model_lib.a nrutil.$(OBJ) locate.$(OBJ) polint.$(OBJ) \
				rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ) zbrent.$(OBJ) \
				hygfx.$(OBJ)  random.$(OBJ) 
locate.$(OBJ)	: locate.f90
	$(FOR) locate.f90 $(FFLAGS)locate.$(OBJ)
polint.$(OBJ)	: polint.f90
	$(FOR) polint.f90 $(FFLAGS)polint.$(OBJ)
nrtype.$(OBJ)	: nrtype.f90
	$(FOR) nrtype.f90 $(FFLAGS)nrtype.$(OBJ)
nr.$(OBJ)	: nr.f90 
	$(FOR) nr.f90 $(FFLAGS)nr.$(OBJ)
nrutil.$(OBJ)	: nrutil.f90
	$(FOR) nrutil.f90 $(FFLAGS)nrutil.$(OBJ)
rkqs.$(OBJ)	: rkqs.f90
	$(FOR) rkqs.f90 $(FFLAGS)rkqs.$(OBJ)	
rkck.$(OBJ)	: rkck.f90
	$(FOR) rkck.f90 $(FFLAGS)rkck.$(OBJ)	
odeint.$(OBJ)	: odeint.f90
	$(FOR) odeint.f90 $(FFLAGS)odeint.$(OBJ)	
zbrent.$(OBJ)	: zbrent.f90
	$(FOR) zbrent.f90 $(FFLAGS2)zbrent.$(OBJ)	
hygfx.$(OBJ) : hygfx.for 
	$(FOR) hygfx.for $(FFLAGS)hygfx.$(OBJ) 
random.$(OBJ) : random.f90 
	$(FOR) random.f90 $(FFLAGS)random.$(OBJ) 
variables.$(OBJ) : variables.f90 nrtype.$(OBJ)
	$(FOR) variables.f90 $(FFLAGS)variables.$(OBJ)
diagnostics.$(OBJ) : diagnostics.f90 mpi_module.$(OBJ) nr.$(OBJ) nrtype.$(OBJ)
	$(FOR) diagnostics.f90 $(FFLAGS)diagnostics.$(OBJ)
initialisation.$(OBJ) : initialisation.f90 random.$(OBJ) nr.$(OBJ) nrtype.$(OBJ) pamm_code
	$(FOR) initialisation.f90 -I ${NETCDFMOD}  $(FFLAGS)initialisation.$(OBJ) \
	    -I$(PAMM_DIR)
driver_code.$(OBJ) : driver_code.f90 nrtype.$(OBJ) dynamics.$(OBJ) \
        sfvt_code sgm_code pamm_code rtm_code lsm_code diagnostics.$(OBJ)
	$(FOR) driver_code.f90 -I ${NETCDFMOD}  $(FFLAGS)driver_code.$(OBJ) -I$(SFVT_DIR) \
	    -I$(SGM_DIR) -I$(PAMM_DIR) -I$(RTM_DIR) -I$(LSM_DIR)
mpi_module.$(OBJ) : mpi_module.f90 
	$(FOR) mpi_module.f90 $(FFLAGS)mpi_module.$(OBJ)
dynamics.$(OBJ) : dynamics.f90 sgm_code rtm_code 
	$(FOR) dynamics.f90 $(FFLAGSOMP)dynamics.$(OBJ) -I$(SGM_DIR) -I$(RTM_DIR)
main.$(OBJ)   : main.f90 variables.$(OBJ) mpi_module.$(OBJ) initialisation.$(OBJ) \
				 driver_code.$(OBJ) dynamics.$(OBJ) pamm_code rtm_code lsm_code
	$(FOR)  main.f90 -I ${NETCDFMOD} $(FFLAGS)main.$(OBJ) -I$(PAMM_DIR) -I$(RTM_DIR) \
	    -I$(LSM_DIR) -I$(RTM_DIR)/pts

sfvt_code:
	$(MAKE) -C $(SFVT_DIR)

sgm_code:
	$(MAKE) -C $(SGM_DIR)

pamm_code:
	$(MAKE) -C $(PAMM_DIR) MPI_PAMM=1

rtm_code:
	$(MAKE) -C $(RTM_DIR)

lsm_code:
	$(MAKE) -C $(LSM_DIR)

clean: 
	rm *.exe  *.o *.mod *~ \
	model_lib.a

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done
	
	
