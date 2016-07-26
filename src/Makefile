.SUFFIXES: .f .F .o .a  .f90 .F90
.PHONY: fdf
include ../arch.make

EXEC = dmft-ed-zt.x
OBJS = mpi.o timer.o constants.o numeric_utils.o \
	   io_units.o dmft_params.o main.o dmft_grid.o dmft_lattice.o \
	   dmft.o ed_basis.o dmft_green.o ground_state.o \
	   impurity_hamiltonian.o lanczos_solver.o arpack_solver.o \
	   eigpair.o debug.o ed_operator.o ed_projection.o

default: all
##################
### FDF Module ###
##################
FDF=fdf/libfdf.a
LIBS+=$(FDF)
FFLAGS+=-I./fdf

####################
### Dependencies ###
####################
$(OBJS): $(FDF)
main.o: mpi.o timer.o io_units.o constants.o \
		numeric_utils.o dmft.o
dmft.o: dmft_lattice.o dmft_grid.o dmft_green.o eigpair.o ground_state.o \
		impurity_hamiltonian.o ed_projection.o
ground_state.o: debug.o lanczos_solver.o arpack_solver.o
arpack_solver.o: impurity_hamiltonian.o
lanczos_solver.o: impurity_hamiltonian.o
impurity_hamiltonian.o: ed_basis.o
dmft_lattice.o: dmft_grid.o
dmft_green.o: eigpair.o ed_operator.o lanczos_solver.o
ed_operator.o: ed_basis.o

$(FDF): 
	cd fdf ; $(MAKE);

all: $(OBJS)
	$(FC) -o $(EXEC) $(LDFLAGS) $(OBJS) $(LIBS) 

clean:
	rm -f *.x *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) clean)

%.o:%.mod
