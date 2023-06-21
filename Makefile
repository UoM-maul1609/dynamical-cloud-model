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

all: main.exe main_ser_1d.exe main_ser_2d.exe

main.exe	:  main.$(OBJ) variables.$(OBJ) mpi_module.$(OBJ) \
			 initialisation.$(OBJ) driver_code.$(OBJ) \
			  model_lib.a  
	$(FOR2) $(FFLAGSOMP)main.exe main.$(OBJ) variables.$(OBJ) mpi_module.$(OBJ) \
		 initialisation.$(OBJ) driver_code.$(OBJ) \
		  -lm model_lib.a -I$(OSNF_DIR) \
		 ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)	 
main_ser_1d.exe	:  main_ser_1d.$(OBJ) variables.$(OBJ) mpi_module.$(OBJ) \
			 initialisation.$(OBJ) driver_code_ser.$(OBJ) \
			  model_lib.a  
	$(FOR2) $(FFLAGSOMP)main_ser_1d.exe main_ser_1d.$(OBJ) variables.$(OBJ) mpi_module.$(OBJ) \
		 initialisation.$(OBJ) driver_code_ser.$(OBJ) \
		  -lm model_lib.a -I$(OSNF_DIR) \
		 ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
main_ser_2d.exe	:  main_ser_2d.$(OBJ) variables.$(OBJ) mpi_module.$(OBJ) \
			 initialisation.$(OBJ) driver_code_ser.$(OBJ) \
			  model_lib.a  
	$(FOR2) $(FFLAGSOMP)main_ser_2d.exe main_ser_2d.$(OBJ) variables.$(OBJ) mpi_module.$(OBJ) \
		 initialisation.$(OBJ) driver_code_ser.$(OBJ) \
		  -lm model_lib.a -I$(OSNF_DIR) \
		 ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
model_lib.a	:   advection_1d.$(OBJ) advection_2d.$(OBJ) advection_3d.$(OBJ) osnf_code 
	$(AR) rc model_lib.a advection_1d.$(OBJ) advection_2d.$(OBJ) advection_3d.$(OBJ) \
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
	$(FOR) variables.f90 -I$(OSNF_DIR) $(FFLAGS)variables.$(OBJ)
initialisation.$(OBJ) : initialisation.f90 osnf_code
	$(FOR) initialisation.f90 -I ${NETCDFMOD} -I$(OSNF_DIR) $(FFLAGS)initialisation.$(OBJ)
driver_code.$(OBJ) : driver_code.f90 advection_1d.$(OBJ) \
advection_2d.$(OBJ) advection_3d.$(OBJ) osnf_code
	$(FOR) driver_code.f90 -I ${NETCDFMOD} -I$(OSNF_DIR) $(FFLAGS)driver_code.$(OBJ)
driver_code_ser.$(OBJ) : driver_code_ser.f90 advection_1d.$(OBJ) \
advection_2d.$(OBJ) advection_3d.$(OBJ) osnf_code
	$(FOR) driver_code_ser.f90 -I ${NETCDFMOD} -I$(OSNF_DIR) $(FFLAGS)driver_code_ser.$(OBJ)
mpi_module.$(OBJ) : mpi_module.f90 osnf_code
	$(FOR) mpi_module.f90 -I$(OSNF_DIR)  -cpp -DVAR_TYPE=$(VAR_TYPE)  $(FFLAGS)mpi_module.$(OBJ)
advection_1d.$(OBJ) : advection_1d.F90 osnf_code
	$(FOR) advection_1d.F90 -cpp -I$(OSNF_DIR) $(FFLAGSOMP)advection_1d.$(OBJ)
advection_2d.$(OBJ) : advection_2d.F90 osnf_code
	$(FOR) advection_2d.F90 -cpp -I$(OSNF_DIR) $(FFLAGSOMP)advection_2d.$(OBJ)
advection_3d.$(OBJ) : advection_3d.F90 osnf_code
	$(FOR) advection_3d.F90 -I$(OSNF_DIR) $(FFLAGSOMP)advection_3d.$(OBJ)
main.$(OBJ)   : main.f90 variables.$(OBJ) mpi_module.$(OBJ) initialisation.$(OBJ) \
				 driver_code.$(OBJ) advection_1d.$(OBJ) advection_2d.$(OBJ) \
				 advection_3d.$(OBJ) 
	$(FOR)  main.f90 -I ${NETCDFMOD} $(FFLAGS)main.$(OBJ) 
main_ser_1d.$(OBJ)   : main_ser_1d.f90 variables.$(OBJ) mpi_module.$(OBJ) initialisation.$(OBJ) \
				 driver_code_ser.$(OBJ) advection_1d.$(OBJ) advection_2d.$(OBJ) \
				 advection_3d.$(OBJ) 
	$(FOR)  main_ser_1d.f90 -I ${NETCDFMOD} $(FFLAGS)main_ser_1d.$(OBJ) 
main_ser_2d.$(OBJ)   : main_ser_2d.f90 variables.$(OBJ) mpi_module.$(OBJ) initialisation.$(OBJ) \
				 driver_code_ser.$(OBJ) advection_1d.$(OBJ) advection_2d.$(OBJ) \
				 advection_3d.$(OBJ) 
	$(FOR)  main_ser_2d.f90 -I ${NETCDFMOD} $(FFLAGS)main_ser_2d.$(OBJ) 

osnf_code:
	$(MAKE) -C $(OSNF_DIR)

clean: 
	rm *.exe  *.o *.mod *~ \
	model_lib.a

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done
	

