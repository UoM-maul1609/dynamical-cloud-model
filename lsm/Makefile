OSNF_DIR = osnf

.PHONY: osnf cleanall
CLEANDIRS = $(OSNF_DIR) ./



DEBUG = -fbounds-check -g 
OPT    =-O3

# these three lines should be edited for your system. On systems 
# that do not have separate fortran and c libraries, set NETCDF_FOR and NETCDF_C
# to the same, and set NETCDF_LIB to -lnetcdf (i.e. without the extra f)
#NETCDF_FOR=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.4-mac/
#NETCDF_C=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.1.1-mac/
NETCDF_LIB=-lnetcdff 

NETCDFLIB=-L ${NETCDF_FOR}/lib/ \
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


main.exe	:  main.$(OBJ) variables.$(OBJ)  mpi_module.$(OBJ) \
			 initialisation.$(OBJ) driver_code.$(OBJ) osnf_code lsm_lib.a 
	$(FOR2) $(FFLAGSOMP)main.exe main.$(OBJ) variables.$(OBJ)  land_surface.$(OBJ) \
			 mpi_module.$(OBJ) \
			 initialisation.$(OBJ) driver_code.$(OBJ) \
			 -lm lsm_lib.a  -I$(OSNF_DIR) \
		 ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
lsm_lib.a	:   land_surface.$(OBJ) osnf_code 
	$(AR) rc lsm_lib.a land_surface.$(OBJ) \
				$(OSNF_DIR)/numerics.$(OBJ) $(OSNF_DIR)/zeroin.$(OBJ) $(OSNF_DIR)/sfmin.$(OBJ) \
                $(OSNF_DIR)/fmin.$(OBJ) $(OSNF_DIR)/r1mach.$(OBJ) \
                $(OSNF_DIR)/d1mach.$(OBJ) $(OSNF_DIR)/dfsid1.$(OBJ) \
                $(OSNF_DIR)/poly_int.$(OBJ) $(OSNF_DIR)/find_pos.$(OBJ) \
                $(OSNF_DIR)/svode.$(OBJ) \
                $(OSNF_DIR)/slinpk.$(OBJ) $(OSNF_DIR)/vode.$(OBJ) \
                $(OSNF_DIR)/dlinpk.$(OBJ) $(OSNF_DIR)/vode_integrate.$(OBJ) \
                $(OSNF_DIR)/erfinv.$(OBJ) $(OSNF_DIR)/tridiagonal.$(OBJ) \
                $(OSNF_DIR)/hygfx.$(OBJ) $(OSNF_DIR)/random.$(OBJ)				
variables.$(OBJ) : variables.f90 osnf_code
	$(FOR) variables.f90 $(FFLAGS)variables.$(OBJ) -I$(OSNF_DIR)
initialisation.$(OBJ) : initialisation.f90 osnf_code
	$(FOR) initialisation.f90 -I ${NETCDFMOD}  $(FFLAGS)initialisation.$(OBJ) -I$(OSNF_DIR)
driver_code.$(OBJ) : driver_code.f90 land_surface.$(OBJ) osnf_code
	$(FOR) driver_code.f90 -I ${NETCDFMOD}  -I$(OSNF_DIR) $(FFLAGS)driver_code.$(OBJ)
mpi_module.$(OBJ) : mpi_module.f90 
	$(FOR) mpi_module.f90  -cpp -DVAR_TYPE=$(VAR_TYPE) $(FFLAGS)mpi_module.$(OBJ) -I$(OSNF_DIR)
land_surface.$(OBJ) : land_surface.f90 mpi_module.$(OBJ) \
					initialisation.$(OBJ) osnf_code
	$(FOR) land_surface.f90  -I$(OSNF_DIR) $(FFLAGSOMP)land_surface.$(OBJ)
main.$(OBJ)   : main.f90 variables.$(OBJ)  osnf_code mpi_module.$(OBJ) \
			initialisation.$(OBJ) land_surface.$(OBJ) driver_code.$(OBJ) 
	$(FOR)  main.f90 -I ${NETCDFMOD} -I$(OSNF_DIR) $(FFLAGS)main.$(OBJ) 

osnf_code:
	$(MAKE) -C $(OSNF_DIR)

clean: 
	rm *.exe  *.o *.mod *~ \
	lsm_lib.a

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done
	
