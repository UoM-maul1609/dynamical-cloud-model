BAM_DIR = bam
SFVT_DIR = sfvt

.PHONY: bam_code cleanall
.PHONY: sfvt_code cleanall
CLEANDIRS = $(BAM_DIR) $(SFVT_DIR) ./


DEBUG = -fbounds-check -g
MPI_PAMM   = 0
OPT    =-O3


FOR = mpif90 -c  
FOR2 = mpif90  

AR = ar 
RANLIB = ranlib 
OBJ = o

FFLAGS = $(OPT)  $(DEBUG) -w -fallow-argument-mismatch -o 
#FFLAGSOMP = -fopenmp-simd $(FFLAGS)
FFLAGS2 =  $(DEBUG) -w -fallow-argument-mismatch -O3 -o 


pmicro_lib.a	:   nrtype.$(OBJ) nr.$(OBJ) nrutil.$(OBJ) locate.$(OBJ) polint.$(OBJ) \
				rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ) zbrent.$(OBJ) dfridr.$(OBJ) \
				hygfx.$(OBJ) microphysics.$(OBJ) advection_1d.$(OBJ) erfinv.$(OBJ) \
				bam_code sfvt_code
	$(AR) rc pmicro_lib.a nrutil.$(OBJ) locate.$(OBJ) polint.$(OBJ) \
				rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ) zbrent.$(OBJ) dfridr.$(OBJ) \
				hygfx.$(OBJ) microphysics.$(OBJ) advection_1d.$(OBJ) erfinv.$(OBJ) \
				$(BAM_DIR)/nr_code.$(OBJ) $(BAM_DIR)/bulk_activation_module.$(OBJ) \
				$(BAM_DIR)/random.$(OBJ) \
				$(SFVT_DIR)/nrutil.$(OBJ) $(SFVT_DIR)/locate.$(OBJ) $(SFVT_DIR)/polint.$(OBJ) \
				$(SFVT_DIR)/rkqs.$(OBJ) $(SFVT_DIR)/rkck.$(OBJ) $(SFVT_DIR)/odeint.$(OBJ) \
				$(SFVT_DIR)/zbrent.$(OBJ) \
				$(SFVT_DIR)/hygfx.$(OBJ) $(SFVT_DIR)/random.$(OBJ) \
				$(SFVT_DIR)/advection_1d.$(OBJ) $(SFVT_DIR)/advection_2d.$(OBJ) \
				$(SFVT_DIR)/advection_3d.$(OBJ)
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
dfridr.$(OBJ)	: dfridr.f90
	$(FOR) dfridr.f90 $(FFLAGS)dfridr.$(OBJ)	
odeint.$(OBJ)	: odeint.f90
	$(FOR) odeint.f90 $(FFLAGS)odeint.$(OBJ)	
zbrent.$(OBJ)	: zbrent.f90
	$(FOR) zbrent.f90 $(FFLAGS2)zbrent.$(OBJ)	
erfinv.$(OBJ) : erfinv.f90
	$(FOR) erfinv.f90 $(FFLAGS)erfinv.$(OBJ) 
advection_1d.$(OBJ) : advection_1d.f90 
	$(FOR) advection_1d.f90 $(FFLAGS)advection_1d.$(OBJ)
microphysics.$(OBJ) : microphysics.f90 advection_1d.$(OBJ) dfridr.$(OBJ) erfinv.$(OBJ) \
                bam_code sfvt_code
	$(FOR) microphysics.f90 -cpp -DMPI_PAMM=$(MPI_PAMM) $(FFLAGS)microphysics.$(OBJ) -I$(BAM_DIR) -I$(SFVT_DIR)
hygfx.$(OBJ) : hygfx.for 
	$(FOR) hygfx.for $(FFLAGS)hygfx.$(OBJ) 

bam_code:
	$(MAKE) -C $(BAM_DIR)
	
sfvt_code:
	$(MAKE) -C $(SFVT_DIR)

clean: 
	rm *.exe  *.o *.mod *~ \
	pmicro_lib.a

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done
	
