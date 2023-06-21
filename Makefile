SFVT_DIR = sfvt
SGM_DIR = sgm
PAMM_DIR = pamm
RTM_DIR = rtm
LSM_DIR = lsm
OSNF_DIR = $(SFVT_DIR)/osnf

.PHONY: sfvt_code sgm_code pamm_code rtm_code lsm_code cleanall
CLEANDIRS = $(SFVT_DIR) $(SGM_DIR) $(PAMM_DIR) $(PAMM_DIR)/sfvt $(PAMM_DIR)/bam \
    $(PAMM_DIR)/sfvt/osnf $(RTM_DIR) $(RTM_DIR)/pts $(LSM_DIR) \
    $(SFVT_DIR)/osnf $(SGM_DIR)/osnf $(PAMM_DIR)/bam/osnf  \
    $(RTM_DIR)/pts/osnf $(LSM_DIR)/osnf ./


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
VAR_TYPE = 1 # 0 single, 1 double


main.exe	:  main.$(OBJ) variables.$(OBJ) mpi_module.$(OBJ) \
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
model_lib.a	:   sfvt_code
	$(AR) rc model_lib.a $(OSNF_DIR)/numerics.$(OBJ) $(OSNF_DIR)/zeroin.$(OBJ) $(OSNF_DIR)/sfmin.$(OBJ) \
				$(OSNF_DIR)/fmin.$(OBJ) $(OSNF_DIR)/r1mach.$(OBJ) \
                $(OSNF_DIR)/d1mach.$(OBJ) $(OSNF_DIR)/dfsid1.$(OBJ) \
                $(OSNF_DIR)/poly_int.$(OBJ) $(OSNF_DIR)/find_pos.$(OBJ) \
                $(OSNF_DIR)/svode.$(OBJ) \
                $(OSNF_DIR)/slinpk.$(OBJ) $(OSNF_DIR)/vode.$(OBJ) \
                $(OSNF_DIR)/dlinpk.$(OBJ) $(OSNF_DIR)/vode_integrate.$(OBJ) \
                $(OSNF_DIR)/erfinv.$(OBJ) $(OSNF_DIR)/tridiagonal.$(OBJ) \
                $(OSNF_DIR)/hygfx.$(OBJ) $(OSNF_DIR)/random.$(OBJ)	
variables.$(OBJ) : variables.f90 sfvt_code
	$(FOR) variables.f90 -I$(OSNF_DIR) $(FFLAGS)variables.$(OBJ)
diagnostics.$(OBJ) : diagnostics.f90 mpi_module.$(OBJ) sfvt_code
	$(FOR) diagnostics.f90 -cpp -DVAR_TYPE=$(VAR_TYPE) -I$(OSNF_DIR) $(FFLAGS)diagnostics.$(OBJ)
initialisation.$(OBJ) : initialisation.f90 sfvt_code pamm_code
	$(FOR) initialisation.f90 -I${NETCDFMOD}  $(FFLAGS)initialisation.$(OBJ) \
	    -I$(PAMM_DIR) -I$(OSNF_DIR) 
driver_code.$(OBJ) : driver_code.f90 dynamics.$(OBJ) \
        sfvt_code sgm_code pamm_code rtm_code lsm_code diagnostics.$(OBJ)
	$(FOR) driver_code.f90 -I${NETCDFMOD}  $(FFLAGS)driver_code.$(OBJ) -I$(SFVT_DIR) \
	    -I$(SGM_DIR) -I$(PAMM_DIR) -I$(RTM_DIR) -I$(LSM_DIR) -I$(OSNF_DIR) 
mpi_module.$(OBJ) : mpi_module.f90 sfvt_code
	$(FOR) mpi_module.f90 -cpp -DVAR_TYPE=$(VAR_TYPE) -I$(OSNF_DIR) $(FFLAGS)mpi_module.$(OBJ)
dynamics.$(OBJ) : dynamics.f90 sgm_code rtm_code sfvt_code
	$(FOR) dynamics.f90 -cpp -DVAR_TYPE=$(VAR_TYPE) \
		$(FFLAGSOMP)dynamics.$(OBJ) -I$(SGM_DIR) -I$(RTM_DIR) -I$(OSNF_DIR) 
main.$(OBJ)   : main.f90 variables.$(OBJ) mpi_module.$(OBJ) initialisation.$(OBJ) \
				 driver_code.$(OBJ) dynamics.$(OBJ) pamm_code rtm_code lsm_code sfvt_code
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
	
	
