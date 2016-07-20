.SUFFIXES: .f .F .o .a  .f90 .F90
.PHONY: fdf
include ../arch.make

EXEC = dmft-ed-zt.x
OBJS = mpi.o timer.o constants.o numeric_utils.o alloc.o\
	   io_units.o dmft_params.o main.o dmft_grid.o dmft_lattice.o

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
main.o: mpi.o timer.o io_units.o alloc.o constants.o\
		numeric_utils.o dmft_grid.o dmft_lattice.o

dmft_lattic.o: dmft_grid.o

$(FDF): 
	cd fdf ; $(MAKE);

all: $(OBJS)
	$(FC) -o $(EXEC) $(LDFLAGS) $(OBJS) $(LIBS) 

clean:
	rm -f *.x *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) clean)

%.o:%.mod
