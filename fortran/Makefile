# vp4d_remapper
# vm4d_remapper
# vm4d_spectral
# vp4d_polar
# vm4d_spectral_charge
# vm2d_spectral_charge
# vp4d_sequential

SIM_NAME = vp4d_sequential

SLL_DIR ?= "/opt/selalib"
FC = h5pfc
FFLAGS = -w -ffree-line-length-none -fall-intrinsics -O3 -fPIC -march=native -I${SLL_DIR}/include/selalib
FLIBS = -L${SLL_DIR}/lib -lselalib -lfftw3 -ldfftpack

${SIM_NAME}: sll_vlasov4d_spectral.o sll_vlasov4d_spectral_charge.o sll_vlasov4d_maxwell.o sll_vlasov4d_poisson.o sll_vlasov4d_polar.o sll_vlasov4d_base.o init_functions.o ${SIM_NAME}.F90
	${FC} ${FFLAGS} -o ${SIM_NAME} $^ -I${SLL_DIR}/include/selalib ${FLIBS}

run: build
	mpirun -np 8 ${SIM_NAME}

clean: 
	rm -f *.o ${SIM_NAME} *.mod

selalib:
	git clone git@gitlab.inria.fr:mingus/selalib.git

sll_build: selalib
	mkdir -p $@
	cd $@; cmake ../selalib -DHDF5_PARALLEL_ENABLED=ON \
                  -DCMAKE_INSTALL_PREFIX=${SLL_DIR} \
                  -DCMAKE_BUILD_TYPE=Release \
 	          -DBUILD_TESTING=OFF -DBUILD_SIMULATIONS=OFF; make install

.SUFFIXES: $(SUFFIXES) .F90

.F90.o:
	$(FC) $(FFLAGS) -c $<

.mod.o:
	$(FC) $(FFLAGS) -c $*.F90
