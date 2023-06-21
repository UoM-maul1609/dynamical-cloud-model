OSNF_DIR = osnf

.PHONY: osnf_code cleanall
CLEANDIRS = $(OSNF_DIR) ./


DEBUG = -fbounds-check -g
MPI    =#-DMPI1
OPT    =-O3


FOR = gfortran -c  
FOR2 = gfortran  

AR = ar 
RANLIB = ranlib 
OBJ = o
FFLAGS = $(OPT)  $(DEBUG)  -o 
FFLAGS2 =  $(DEBUG) $(OPT) -o 
VAR_TYPE = 1 # 0 single, 1 double


main.exe	:  bam_lib.a  main.$(OBJ) bulk_activation_module.$(OBJ)
	$(FOR2) $(FFLAGS2)main.exe main.$(OBJ) -lm bam_lib.a $(DEBUG)
bam_lib.a	:   bulk_activation_module.$(OBJ) 
	$(AR) rc bam_lib.a bulk_activation_module.$(OBJ)  \
				$(OSNF_DIR)/numerics.$(OBJ) $(OSNF_DIR)/zeroin.$(OBJ) $(OSNF_DIR)/sfmin.$(OBJ) \
                $(OSNF_DIR)/fmin.$(OBJ) $(OSNF_DIR)/r1mach.$(OBJ) \
                $(OSNF_DIR)/d1mach.$(OBJ) $(OSNF_DIR)/dfsid1.$(OBJ) \
                $(OSNF_DIR)/poly_int.$(OBJ) $(OSNF_DIR)/find_pos.$(OBJ) \
                $(OSNF_DIR)/svode.$(OBJ) \
                $(OSNF_DIR)/slinpk.$(OBJ) $(OSNF_DIR)/vode.$(OBJ) \
                $(OSNF_DIR)/dlinpk.$(OBJ) $(OSNF_DIR)/vode_integrate.$(OBJ) \
                $(OSNF_DIR)/erfinv.$(OBJ) $(OSNF_DIR)/tridiagonal.$(OBJ) \
                $(OSNF_DIR)/hygfx.$(OBJ) $(OSNF_DIR)/random.$(OBJ)				
	
main.$(OBJ)	: main.f90 osnf_code
	$(FOR) main.f90 $(FFLAGS)main.$(OBJ) -I$(OSNF_DIR)
bulk_activation_module.$(OBJ)	: bulk_activation_module.f90 osnf_code
	$(FOR) bulk_activation_module.f90 $(FFLAGS)bulk_activation_module.$(OBJ) -I$(OSNF_DIR)
osnf_code:
	$(MAKE) -C $(OSNF_DIR)

clean: 
	rm *.exe  *.o *.mod *~ \
	bam_lib.a

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done
