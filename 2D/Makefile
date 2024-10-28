FC = mpif90
FFLAGS=-O3 -Wunused-variable -g -fbacktrace 
mcmodel=large

SRCS := src/main.f90 \
	src_gradient/checkpointing_double.f90 \
	src_gradient/compute_adjoint.f90 \
	src_gradient/derivative_first_line_moinsdegen.f90 \
	src_gradient/compute_kernel_interpolation.f90 \
	src/computePML.f90 \
	src/discretisation.f90 \
	src_forward/forwardproblem.f90 \
	src/initialisation.f90 \
	src/mpi_comm.f90 \
	src/outputs_lib.f90  \
	src/reset.f90 \
	src_inversion/misfit.f90 \
	src_inversion/inversion.f90 \


OBJS := ${SRCS:.f90=.o} 



xinfrasound: parameters.o  $(OBJS)
	$(FC) $(FFLAGS) -mcmodel=${mcmodel} -o xinfrasound parameters.o $(OBJS)
	
%.o: %.f90
	$(FC) $(FFLAGS) -mcmodel=${mcmodel} -o $@ -c $<
	
parameters.o: parameters.f90
	$(FC) $(FFLAGS) -c parameters.f90

clean:
	@rm -f *.mod *.o xinfrasound src_gradient/*.o src/*.o src_forward/*.o src_inversion/*.o
