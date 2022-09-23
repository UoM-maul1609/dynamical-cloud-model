.PHONY: cleanall
CLEANDIRS =  ./


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


main.exe	:  main.$(OBJ) variables.$(OBJ)  mpi_module.$(OBJ) \
			 initialisation.$(OBJ) driver_code.$(OBJ) lsm_lib.a 
	$(FOR2) $(FFLAGSOMP)main.exe main.$(OBJ) variables.$(OBJ)  land_surface.$(OBJ) \
			 mpi_module.$(OBJ) \
			 initialisation.$(OBJ) driver_code.$(OBJ) \
			 -lm lsm_lib.a \
		 ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
lsm_lib.a	:   land_surface.$(OBJ) \
                nrtype.$(OBJ) nr.$(OBJ) nrutil.$(OBJ) locate.$(OBJ) polint.$(OBJ) \
				trapzd.$(OBJ) qromb.$(OBJ) \
				tridag.$(OBJ) rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ) zbrent.$(OBJ) \
				random.$(OBJ) 
	$(AR) rc lsm_lib.a land_surface.$(OBJ) \
	            nrutil.$(OBJ) locate.$(OBJ) polint.$(OBJ) trapzd.$(OBJ) qromb.$(OBJ) \
	            tridag.$(OBJ) \
				rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ) zbrent.$(OBJ) \
				random.$(OBJ) 
locate.$(OBJ)	: locate.f90
	$(FOR) locate.f90 $(FFLAGS)locate.$(OBJ)
polint.$(OBJ)	: polint.f90 nrutil.$(OBJ)
	$(FOR) polint.f90 $(FFLAGS)polint.$(OBJ)
trapzd.$(OBJ)	: trapzd.f90
	$(FOR) trapzd.f90 $(FFLAGS)trapzd.$(OBJ)
qromb.$(OBJ)	: qromb.f90
	$(FOR) qromb.f90 $(FFLAGS)qromb.$(OBJ)
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
	$(FOR) initialisation.f90 -I ${NETCDFMOD}  $(FFLAGS)initialisation.$(OBJ)
driver_code.$(OBJ) : driver_code.f90 nrtype.$(OBJ) land_surface.$(OBJ) 
	$(FOR) driver_code.f90 -I ${NETCDFMOD}  $(FFLAGS)driver_code.$(OBJ)
mpi_module.$(OBJ) : mpi_module.f90 
	$(FOR) mpi_module.f90 $(FFLAGS)mpi_module.$(OBJ)
land_surface.$(OBJ) : land_surface.f90 nrtype.$(OBJ) nr.$(OBJ) mpi_module.$(OBJ) \
					initialisation.$(OBJ) locate.$(OBJ) polint.$(OBJ) \
					trapzd.$(OBJ) qromb.$(OBJ) 
	$(FOR) land_surface.f90 $(FFLAGSOMP)land_surface.$(OBJ) 
main.$(OBJ)   : main.f90 variables.$(OBJ) nrtype.$(OBJ)  mpi_module.$(OBJ) \
			initialisation.$(OBJ) land_surface.$(OBJ) driver_code.$(OBJ) 
	$(FOR)  main.f90 -I ${NETCDFMOD} $(FFLAGS)main.$(OBJ) 

clean: 
	rm *.exe  *.o *.mod *~ \
	lsm_lib.a

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done
	
