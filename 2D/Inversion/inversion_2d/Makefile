FC = mpif90
FFLAGS=-O3

SRCS := main.f90 \
	misfit.f90 \
	inversion.f90 \
	gradient/checkpointing_double.f90 \
	gradient/compute_adjoint.F90 \
	gradient/compute_kernel_interpolation.F90 \
	gradient/computePML.f90 \
	gradient/discretisation.f90 \
	gradient/forwardproblem.f90 \
	gradient/initialisation.f90 \
	gradient/mpi_comm.f90 \
	gradient/outputs_lib.f90  \
	gradient/reset.f90 


OBJS := ${SRCS:.f90=.o} 


inversion: parameters.o  $(OBJS)
	$(FC) $(FFLAGS) -o inversion parameters.o $(OBJS)
	
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
	
parameters.o: parameters.f90
	$(FC) $(FFLAGS) -c parameters.f90

clean:
	@rm -f *.mod *.o inversion gradient/*.o
