

CC = gcc
CFLAGS +=  -I. -O3 -std=c99 -Wall -Wextra
 
LINKLIB = gfortran -shared
#LIBS =  -L$(MULTINEST) -lnest3 -llapack -lcblas -lf77blas -latlas -lpthread
LIBS =  -L$(MULTINEST) -lnest3 -llapack -lcblas -latlas -lpthread

libcnest.so: $(OBJFILES)
	# 
	# MULTINEST set to ${MULTINEST}. 
	# 
	# If you didn't set this to the directory where MultiNest 
	# is compiled, you will receive see this error later:
	#     OSError: libcnest.so: undefined symbol: __nested_MOD_nestrun
	#
	# Checking if libnest3.so is in $MULTINEST
	test -e ${MULTINEST}/libnest3.so
	# Found.
	#
	$(LINKLIB) $(LIBS) $(CFLAGS) -fPIC -std=c99 -shared cnest.c -o libcnest.so
	# all done.

clean:
	rm -f *.so
