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


main.exe	:  bam_lib.a  main.$(OBJ) bulk_activation_module.$(OBJ) nr_code.$(OBJ)
	$(FOR2) $(FFLAGS2)main.exe main.$(OBJ) -lm bam_lib.a $(DEBUG)
bam_lib.a	:   nr_code.$(OBJ) bulk_activation_module.$(OBJ) random.$(OBJ)
	$(AR) rc bam_lib.a nr_code.$(OBJ) bulk_activation_module.$(OBJ) random.$(OBJ)
main.$(OBJ)	: main.f90
	$(FOR) main.f90 $(FFLAGS)main.$(OBJ) 
nr_code.$(OBJ)	: nr_code.f90
	$(FOR) nr_code.f90 $(FFLAGS)nr_code.$(OBJ)
bulk_activation_module.$(OBJ)	: bulk_activation_module.f90
	$(FOR) bulk_activation_module.f90 $(FFLAGS)bulk_activation_module.$(OBJ)
random.$(OBJ) : random.f90 
	$(FOR) random.f90 $(FFLAGS)random.$(OBJ) 
clean :
	rm *.exe  *.o *.mod *~ \
	bam_lib.a

