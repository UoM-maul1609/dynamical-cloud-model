BAM_DIR = bam
SFVT_DIR = sfvt
OSNF_DIR = $(SFVT_DIR)/osnf

.PHONY: bam_code cleanall
.PHONY: sfvt_code cleanall
CLEANDIRS = $(BAM_DIR) $(SFVT_DIR) $(BAM_DIR)/osnf $(SFVT_DIR)/osnf ./


DEBUG = -fbounds-check -g
MPI_PAMM   = 0
OPT    =-O3


FOR = mpif90 -c  
FOR2 = mpif90  

AR = ar 
RANLIB = ranlib 
OBJ = o

FFLAGS = $(OPT)  $(DEBUG) -w -o 
#FFLAGSOMP = -fopenmp-simd $(FFLAGS)
FFLAGS2 =  $(DEBUG) -w -O3 -o 
VAR_TYPE = 1 # 0 single, 1 double


pmicro_lib.a	:   microphysics.$(OBJ) advection_1d.$(OBJ)  \
				bam_code sfvt_code
	$(AR) rc pmicro_lib.a microphysics.$(OBJ) advection_1d.$(OBJ) \
				$(BAM_DIR)/bulk_activation_module.$(OBJ) \
				$(SFVT_DIR)/advection_1d.$(OBJ) $(SFVT_DIR)/advection_2d.$(OBJ) \
				$(SFVT_DIR)/advection_3d.$(OBJ) \
				$(OSNF_DIR)/numerics.$(OBJ) $(OSNF_DIR)/zeroin.$(OBJ) $(OSNF_DIR)/sfmin.$(OBJ) \
                $(OSNF_DIR)/fmin.$(OBJ) $(OSNF_DIR)/r1mach.$(OBJ) \
                $(OSNF_DIR)/d1mach.$(OBJ) $(OSNF_DIR)/dfsid1.$(OBJ) \
                $(OSNF_DIR)/poly_int.$(OBJ) $(OSNF_DIR)/find_pos.$(OBJ) \
                $(OSNF_DIR)/svode.$(OBJ) \
                $(OSNF_DIR)/slinpk.$(OBJ) $(OSNF_DIR)/vode.$(OBJ) \
                $(OSNF_DIR)/dlinpk.$(OBJ) $(OSNF_DIR)/vode_integrate.$(OBJ) \
                $(OSNF_DIR)/erfinv.$(OBJ) $(OSNF_DIR)/tridiagonal.$(OBJ) \
                $(OSNF_DIR)/hygfx.$(OBJ) $(OSNF_DIR)/random.$(OBJ)				
advection_1d.$(OBJ) : advection_1d.f90 sfvt_code
	$(FOR) advection_1d.f90 -I$(OSNF_DIR) $(FFLAGS)advection_1d.$(OBJ)
microphysics.$(OBJ) : microphysics.f90 advection_1d.$(OBJ)  \
                bam_code sfvt_code 
	$(FOR) microphysics.f90 -cpp -DMPI_PAMM=$(MPI_PAMM) $(FFLAGS)microphysics.$(OBJ) \
	-I$(BAM_DIR) -I$(SFVT_DIR) -I$(OSNF_DIR)

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
	
