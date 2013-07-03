
# Change next the right multinest library name
# (tipical are: nest3, multinest, multinest_mpi, ...)
MNESTLIBNAME := multinest
# MNESTLIBNAME := nest3

CC = gcc
FC = gfortran
FFLAGS +=  -I. -O3 -w -fPIC
CFLAGS +=  -I. -O3 -std=c99 -Wall -Wextra -fPIC

# If the MultiNest function is named 
# differently; then do e.g.
# $ make MULTINEST_CALL=__nested__nestrun

ifneq ($(shell which mpif90),)
ifndef WITHOUT_MPI
# compiling with MPI support. You can disable this with
# $ make WITHOUT_MPI=1
FC = mpif90
CC = mpicc
FFLAGS += -DMPI
CFLAGS += -DMPI
  ifndef MULTINEST_CALL
     MULTINEST_CALL := __nested_MOD_nestrun
     #__nested_nestrun
  endif
endif
else
  $(warning "mpif90 not found, so will compile without MPI!")
endif

ifndef MULTINEST_CALL
  MULTINEST_CALL := __nested_MOD_nestrun
endif

CFLAGS += -DMULTINEST_CALL=${MULTINEST_CALL}

#LINKLIB = ${FC} -shared
LINKLIB = ${CC} -shared
LIBS =  -L$(MULTINEST) -l$(MNESTLIBNAME) -llapack -lpthread
# sometimes, this is needed:
# LIBS +=  -lcblas -latlas
# LIBS +=  -lcblas -lf77blas -latlas

libcnest.so: $(OBJFILES)
	#
	# MULTINEST (the library directory) is set to ${MULTINEST}.
	#
	# MNESTLIBNAME preset to '${MNESTLIBNAME}'. Change this according
	#  to your MultiNest installation.
	#
	# If you didn't set this to the directory where MultiNest 
	# is compiled, you will receive see this error later:
	#     OSError: libcnest.so: undefined symbol: __nested_MOD_nestrun
	#
	# Checking if lib${MNESTLIBNAME}.so is in $$MULTINEST...
	@if [ ! -e ${MULTINEST}/lib${MNESTLIBNAME}.so ]; then echo -e "\nERROR:  lib'$$MNESTLIBNAME'.so not found in '$$MULTINEST'. \nERROR:  Set MULTINEST to the directory where MultiNest is, i.e. where you compiled libmultinest.so.\nERROR:  manual: http://johannesbuchner.github.com/PyMultiNest/install.html\n"; exit 1; fi
	# ...Found!
	#
	$(LINKLIB) $(CFLAGS) cnest.c -o libcnest.so $(LIBS)
	# all done.

clean:
	rm -f *.so

