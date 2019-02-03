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
FFLAGS = $(OPT)  $(DEBUG)  -o 
FFLAGSOMP = -fopenmp-simd $(FFLAGS)
FFLAGS2 =  $(DEBUG) -O3 -o 


main.exe	:  main.$(OBJ) variables.$(OBJ)  parallel_tridag.$(OBJ) mpi_module.$(OBJ) \
			 initialisation.$(OBJ) pts_lib.a 
	$(FOR2) $(FFLAGSOMP)main.exe main.$(OBJ) variables.$(OBJ)  parallel_tridag.$(OBJ) \
			 mpi_module.$(OBJ) \
			 initialisation.$(OBJ) -lm pts_lib.a \
		     $(DEBUG)
# 	$(FOR2) $(FFLAGSOMP)main.exe main.$(OBJ) variables.$(OBJ)  parallel_tridag.$(OBJ) \
# 			 mpi_module.$(OBJ) \
# 			 initialisation.$(OBJ) -lm model_lib.a \
# 		 ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
pts_lib.a	:   nrtype.$(OBJ) nr.$(OBJ) nrutil.$(OBJ) locate.$(OBJ) polint.$(OBJ) \
				tridag.$(OBJ) rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ) zbrent.$(OBJ) \
				random.$(OBJ)
	$(AR) rc pts_lib.a nrutil.$(OBJ) locate.$(OBJ) polint.$(OBJ) tridag.$(OBJ) \
				rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ) zbrent.$(OBJ) \
				random.$(OBJ)
locate.$(OBJ)	: locate.f90
	$(FOR) locate.f90 $(FFLAGS)locate.$(OBJ)
polint.$(OBJ)	: polint.f90
	$(FOR) polint.f90 $(FFLAGS)polint.$(OBJ)
tridag.$(OBJ)	: tridag.f90
	$(FOR) tridag.f90 $(FFLAGS)tridag.$(OBJ)
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
random.$(OBJ) : random.f90 
	$(FOR) random.f90 $(FFLAGS)random.$(OBJ) 
variables.$(OBJ) : variables.f90 nrtype.$(OBJ)
	$(FOR) variables.f90 $(FFLAGS)variables.$(OBJ)
initialisation.$(OBJ) : initialisation.f90 random.$(OBJ) nr.$(OBJ) nrtype.$(OBJ)
	$(FOR) initialisation.f90 $(FFLAGS)initialisation.$(OBJ)
# 	$(FOR) initialisation.f90 -I ${NETCDFMOD}  $(FFLAGS)initialisation.$(OBJ)
mpi_module.$(OBJ) : mpi_module.f90 
	$(FOR) mpi_module.f90 $(FFLAGS)mpi_module.$(OBJ)
parallel_tridag.$(OBJ) : parallel_tridag.f90 nrtype.$(OBJ) nr.$(OBJ) mpi_module.$(OBJ) \
					initialisation.$(OBJ)
	$(FOR) parallel_tridag.f90 $(FFLAGSOMP)parallel_tridag.$(OBJ)
main.$(OBJ)   : main.f90 variables.$(OBJ) nrtype.$(OBJ)  mpi_module.$(OBJ) \
			initialisation.$(OBJ) parallel_tridag.$(OBJ) 
	$(FOR)  main.f90 $(FFLAGS)main.$(OBJ) 
# 	$(FOR)  main.f90 -I ${NETCDFMOD} $(FFLAGS)main.$(OBJ) 

clean :
	rm *.exe  *.o *.mod *~ \
	pts_lib.a

