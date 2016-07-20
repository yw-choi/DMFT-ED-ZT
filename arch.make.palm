
FC      = /opt/mvapich2-x/intel/bin/mpif90
FFLAGS = -traceback -fast -no-ipo -xSSE4.2 -I$(MKLROOT)/include/fftw

# FFLAGS  = -traceback -C -CB -I$(MKLROOT)/include/fftw -debug all
FFTW    = $(MKLROOT)/interfaces/fftw3xf/libfftw3xf_intel.a
ARPACK  = /lily/ywchoe/workspace/DMFT-ED/lib/ARPACK/libarpack_LINUX.a
PARPACK = /lily/ywchoe/workspace/DMFT-ED/lib/ARPACK/parpack_MPI-LINUX.a
MKL     = -mkl=sequential
LIBS    = $(MKL) $(PARPACK) $(ARPACK) $(FFTW)

################################################################
DEFS            = -DMPI 
#
.F.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS)   $<
#
