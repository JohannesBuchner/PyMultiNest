Installing PyMultiNest and PyCuba
=================================================

.. contents::

You need to install the python module and put the libraries it uses
into your library path.

1. Installing the Python Module
---------------------------------

Installing the python module from PyPI is easy:

.. code-block:: bash

	$ pip install pymultinest

* Use the "--user" switch if you only want to install the software locally
* On older systems, you may need to use easy_install instead of "pip install"

To get the latest version, download the source directly from the
`Code repository <https://github.com/JohannesBuchner/PyMultiNest/>`_
and install it locally:

.. code-block:: bash

	$ git clone https://github.com/JohannesBuchner/PyMultiNest/
	$ cd PyMultiNest
	$ python setup.py install


However, installing pymultinest is not enough; MultiNest and Cuba are not included, so you
will see an error like::

	ERROR:   Could not load MultiNest library "libmultinest.so"
	ERROR:   You have to build it first, and point the LD_LIBRARY_PATH environment variable to it!

or::

	OSError: libcuba.so: cannot open shared object file: No such file or directory

You need to install the libraries themselves! The next section explains how.

2. Prerequisites for building the libraries
-------------------------------------------------------

To use PyMultiNest and PyCuba, you need to install the relevant
libraries. If you only want to use one of them, skip the others according.

* Prerequisites: *numpy*, *scipy*, *matplotlib*, c and fortran compilers (e.g. *gcc*, *gfortran*)
* Prerequisites for MultiNest: *cmake*, *blas*, *lapack*, *atlas* (and their development versions)
* Recommended: *git*

On **Ubuntu Linux**, install with:

.. code-block:: bash

	$ sudo apt-get install python-{scipy,numpy,matplotlib,progressbar} ipython libblas{3,-dev} liblapack{3,-dev} libatlas{3-base,-dev} cmake build-essential git gfortran

On **Mac OSX**:

* As in the instructions above, you need cmake (e.g. with "brew install cmake"), a Fortran compiler (e.g. with "brew install gcc") and possibly MPI (e.g. with "brew install open-mpi" and then "pip install mpi4py")
* If you google for "MultiNest Mac OSX" or "PyMultiNest Mac OSX" you will find installation instructions.

2. Building the libraries
--------------------------------------

* If you want to use :doc:`PyMultiNest <pymultinest>`:

  * Get and compile MultiNest (use the cmake version from `https://github.com/JohannesBuchner/MultiNest`). The goal is to create lib/libmultinest.so ::
  	
  	git clone https://github.com/JohannesBuchner/MultiNest
  	cd MultiNest/build
  	cmake ..
  	make
    
  * On e.g. Mac OSX, make sure the correct compilers are used by calling cmake using -DCMAKE_C_COMPILER=/path/to/gcc -DCMAKE_CXX_COMPILER=/path/to/g++
  * Include the lib/ directory in your LD_LIBRARY_PATH
  * More detailed install instructions for MultiNest are available in the `tutorial <http://johannesbuchner.github.io/pymultinest-tutorial/install.html#on-your-own-computer>`_.

* If you want to use :doc:`PyCuba <pycuba>`:
  Get and compile `Cuba <http://www.feynarts.de/cuba/>`_. The goal is to create the libcuba.so file ::
  
	git clone https://github.com/JohannesBuchner/cuba/
	cd cuba
	./configure
	./makesharedlib.sh
  
  Include the containing directory your LD_LIBRARY_PATH.

* To install on Mac OSX, make sure the multinest/cuba libraries are in your library path. For that, it is simplest to copy them, e.g.::

	$ cp -v ~/Downloads/MultiNest/lib/lib* /anaconda3/lib/

A discussion on installing on Mac can be found in `issue 10 <https://github.com/JohannesBuchner/PyMultiNest/issues/10>`_. Compiling with MPI support on Mac is discussed in `issue 45 <https://github.com/JohannesBuchner/PyMultiNest/issues/45>`_

3. Running some code
--------------------------

PyMultiNest and PyCuba have to be able to find the corresponding 
libraries. So put the three directories in the dynamic library load path:

.. code-block:: bash

     $ export LD_LIBRARY_PATH=$HOME/Downloads/MultiNest/lib:$HOME/Downloads/cuba/directory/:$LD_LIBRARY_PATH

* On Mac OSX, do the same for DYLD_LIBRARY_PATH.
* Replace the above with your actual path.
* Consider putting this line into your shell startup script (e.g. ~/.bashrc).

Test importing the libraries:

.. code-block:: bash

     $ python -c 'import pymultinest'
     $ python -c 'import pycuba'

Try out the demo programs distributed in the package:

.. code-block:: bash

     $ python $OLDPWD/pymultinest_demo.py
        ....
	Acceptance Rate:                        0.690765
	Replacements:                               3650
	Total Samples:                              5284
	Nested Sampling ln(Z):                235.562844
	Importance Nested Sampling ln(Z):     236.164929 +/-  0.147246
	Acceptance Rate:                        0.690809
	Replacements:                               3653
	Total Samples:                              5288
	Nested Sampling ln(Z):                235.565469
	Importance Nested Sampling ln(Z):     236.165091 +/-  0.147221
	 ln(ev)=   235.91594564793959      +/-  0.12311459261215110     
	 Total Likelihood Evaluations:         5288
	 Sampling finished. Exiting MultiNest
	  analysing data from chains/3-.txt

	evidence: 235.9 +- 0.1

	parameter values:
		      x : 15.968 +- 8.548
		      y : 15.165 +- 9.195
     
     $ python $OLDPWD/pymultinest_solver_demo.py
     
     $ python $OLDPWD/pycuba_demo.py

Congratulations! You are now ready to run your own code. Copy the demo files as starting points, play with the functions and analysis, and integrate it to your own code. The documentation should help you:

* Continue with :doc:`PyMultiNest documentation <pymultinest>`
* Continue with :doc:`PyCuba documentation <pycuba>`

Getting help
----------------------------

Try searching the error message. Search through the `existing questions <https://github.com/JohannesBuchner/PyMultiNest/issues?q=>`_.

Generating the documentation
----------------------------

Go in the doc directory and run make:

.. code-block:: bash

     $ cd doc && make html

Point your web browser to _build/html/index.html in doc.

