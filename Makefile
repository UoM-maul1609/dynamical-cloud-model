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
FFLAGS = $(OPT)  $(DEBUG)  -o 
FFLAGSOMP = -fopenmp-simd $(FFLAGS)
FFLAGS2 =  $(DEBUG) -O3 -o 

all: main.exe main_ser_1d.exe main_ser_2d.exe

main.exe	:  main.$(OBJ) variables.$(OBJ) nrtype.$(OBJ) mpi_module.$(OBJ) \
			 initialisation.$(OBJ) driver_code.$(OBJ) \
			  model_lib.a  
	$(FOR2) $(FFLAGSOMP)main.exe main.$(OBJ) variables.$(OBJ) mpi_module.$(OBJ) \
		 initialisation.$(OBJ) driver_code.$(OBJ) \
		  -lm model_lib.a \
		 ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)	 
main_ser_1d.exe	:  main_ser_1d.$(OBJ) variables.$(OBJ) nrtype.$(OBJ) mpi_module.$(OBJ) \
			 initialisation.$(OBJ) driver_code_ser.$(OBJ) \
			  model_lib.a  
	$(FOR2) $(FFLAGSOMP)main_ser_1d.exe main_ser_1d.$(OBJ) variables.$(OBJ) mpi_module.$(OBJ) \
		 initialisation.$(OBJ) driver_code_ser.$(OBJ) \
		  -lm model_lib.a \
		 ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
main_ser_2d.exe	:  main_ser_2d.$(OBJ) variables.$(OBJ) nrtype.$(OBJ) mpi_module.$(OBJ) \
			 initialisation.$(OBJ) driver_code_ser.$(OBJ) \
			  model_lib.a  
	$(FOR2) $(FFLAGSOMP)main_ser_2d.exe main_ser_2d.$(OBJ) variables.$(OBJ) mpi_module.$(OBJ) \
		 initialisation.$(OBJ) driver_code_ser.$(OBJ) \
		  -lm model_lib.a \
		 ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
model_lib.a	:   nrtype.$(OBJ) nr.$(OBJ) nrutil.$(OBJ) locate.$(OBJ) polint.$(OBJ) \
				rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ) zbrent.$(OBJ) \
				hygfx.$(OBJ)  random.$(OBJ) advection_1d.$(OBJ) advection_2d.$(OBJ) \
				advection_3d.$(OBJ)
	$(AR) rc model_lib.a nrutil.$(OBJ) locate.$(OBJ) polint.$(OBJ) \
				rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ) zbrent.$(OBJ) \
				hygfx.$(OBJ)  random.$(OBJ) advection_1d.$(OBJ) advection_2d.$(OBJ) \
				advection_3d.$(OBJ)
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
initialisation.$(OBJ) : initialisation.f90 random.$(OBJ) nr.$(OBJ) nrtype.$(OBJ)
	$(FOR) initialisation.f90 -I ${NETCDFMOD}  $(FFLAGS)initialisation.$(OBJ)
driver_code.$(OBJ) : driver_code.f90 nrtype.$(OBJ) advection_1d.$(OBJ) \
advection_2d.$(OBJ) advection_3d.$(OBJ) 
	$(FOR) driver_code.f90 -I ${NETCDFMOD}  $(FFLAGS)driver_code.$(OBJ)
driver_code_ser.$(OBJ) : driver_code_ser.f90 nrtype.$(OBJ) advection_1d.$(OBJ) \
advection_2d.$(OBJ) advection_3d.$(OBJ) 
	$(FOR) driver_code_ser.f90 -I ${NETCDFMOD}  $(FFLAGS)driver_code_ser.$(OBJ)
mpi_module.$(OBJ) : mpi_module.f90 
	$(FOR) mpi_module.f90 $(FFLAGS)mpi_module.$(OBJ)
advection_1d.$(OBJ) : advection_1d.F90 
	$(FOR) advection_1d.F90 -cpp $(FFLAGSOMP)advection_1d.$(OBJ)
advection_2d.$(OBJ) : advection_2d.F90 
	$(FOR) advection_2d.F90 -cpp $(FFLAGSOMP)advection_2d.$(OBJ)
advection_3d.$(OBJ) : advection_3d.F90 
	$(FOR) advection_3d.F90 $(FFLAGSOMP)advection_3d.$(OBJ)
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

clean :
	rm *.exe  *.o *.mod *~ \
	model_lib.a

