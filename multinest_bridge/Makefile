

CC = gcc
FC = gfortran
FFLAGS +=  -I. -O3 -w -fPIC
CFLAGS +=  -I. -O3 -std=c99 -Wall -Wextra -fPIC
ifndef MULTINEST_CALL
# sometimes it is __nested__nestrun
# then you can run 
# $ MULTINEST_CALL=__nested__nestrun make
# or 
# $ make MULTINEST_CALL=__nested__nestrun
#
MULTINEST_CALL := __nested_MOD_nestrun
endif

ifneq ($(shell which mpif90),)
FC = mpif90
CC = mpicc
FFLAGS += -DMPI
CFLAGS += -DMPI
else
$(warning "mpif90 not found, so will compile without MPI!")
endif

CFLAGS += -DMULTINEST_CALL=${MULTINEST_CALL}

LINKLIB = ${FC} -shared
LIBS =  -L$(MULTINEST) -lnest3 -llapack -lpthread
# sometimes, this is needed:
# LIBS +=  -lcblas -lf77blas -latlas

libcnest.so: $(OBJFILES)
	# 
	# MULTINEST set to ${MULTINEST}. 
	# 
	# If you didn't set this to the directory where MultiNest 
	# is compiled, you will receive see this error later:
	#     OSError: libcnest.so: undefined symbol: __nested_MOD_nestrun
	#
	# Checking if libnest3.so is in $$MULTINEST
	test -e ${MULTINEST}/libnest3.so
	# Found.
	#
	$(LINKLIB) $(LIBS) $(CFLAGS) cnest.c -o libcnest.so
	# all done.

clean:
	rm -f *.so
