Installing PyMultiNest and PyAPEMoST
=================================================

.. contents::

You need to install the python module and put the libraries it uses
into your library path.

1. Installing the Python Module
---------------------------------

Installing the python module from PyPI is easy:

.. code-block:: bash

	$ pip install pymultinest

On older systems, you may need to use easy_install instead of "pip install"
Use the "--user" switch if you only want to install the software locally

To get the latest version, download the source directly from the
`Code repository <https://github.com/JohannesBuchner/PyMultiNest/>`_
and install it locally:

.. code-block:: bash

	$ git clone https://github.com/JohannesBuchner/PyMultiNest/
	$ cd PyMultiNest
	$ python setup.py install

2. Installing the libraries on Linux
--------------------------------------

To use PyMultiNest, PyAPEMoST, and PyCuba, you need to install the relevant
libraries. If you only want to use one of them, skip the others according.

* Prerequisites: *numpy*, *scipy*, *matplotlib*, c and fortran compilers (e.g. *gcc*, *gfortran*)
* Prerequisites for MultiNest: *cmake*, *blas*, *lapack*, *atlas* (and their development versions)
* Recommended: *git*

On Ubuntu, install with:

.. code-block:: bash

	$ sudo apt-get install python-{scipy,numpy,matplotlib,progressbar} ipython libblas{3,-dev} liblapack{3,-dev} libatlas{3-base,-dev} cmake build-essential git gfortran

* If you want to use :doc:`PyMultiNest <pymultinest>`:
  Get and compile MultiNest (use the cmake version from `here <https://github.com/JohannesBuchner/MultiNest>`_). You will have a lib/libmultinest.so. 
  Include the lib/ directory in your LD_LIBRARY_PATH.
  More detailed install instructions for MultiNest are available in the `tutorial <http://johannesbuchner.github.io/pymultinest-tutorial/install.html#on-your-own-computer>`_.

* If you want to use :doc:`PyCuba <pycuba>`:
  Get and compile `Cuba <http://www.feynarts.de/cuba/>`_. You will have a libcuba.so file.
  Include the containing directory your LD_LIBRARY_PATH.

* If you want to use :doc:`PyAPEMoST <pyapemost>`:
  Get and compile `APEMoST <http://apemost.sourceforge.net/doc/>`_. You will have a libapemost.so file.
  Include the containing directory your LD_LIBRARY_PATH.

Running some code
--------------------------

PyMultiNest, PyAPEMoST and PyCuba have to be able to find the corresponding 
libraries. So put the three directories in the dynamic library load path:

.. code-block:: bash

     $ export LD_LIBRARY_PATH=/my/multinest/directory/lib:/my/cuba/directory/:/my/apemost/directory/:$LD_LIBRARY_PATH

Consider putting this line into your shell startup script (e.g. ~/.bashrc).

Test importing the libraries:

.. code-block:: bash

     $ python -c 'import pymultinest'
     $ python -c 'import pyapemost'

Try out the demo programs distributed in the package. They produce a lot of output, so lets create a temporary directory:

.. code-block:: bash

     $ mkdir /tmp/demo_output
     $ cd /tmp/demo_output
     
     $ python $OLDPWD/pymultinest_demo.py
     
     $ python $OLDPWD/pyapemost_demo.py

Congratulations! You are now ready to run your own code. Copy the demo files as starting points, play with the functions and analysis, and integrate it to your own code. The documentation should help you.

Getting help
----------------------------
Try searching the error message. Look at the `Open issues <https://github.com/JohannesBuchner/PyMultiNest/issues?state=open>`_.

Generating the documentation
----------------------------

Go in the doc directory and run make:

.. code-block:: bash

     $ cd doc && make html

Point your web browser to _build/html/index.html in doc.


Install on Mac
-----------------------------
Follow the above instructions. Use DYLD_LIBRARY_PATH instead of LD_LIBRARY_PATH.
rename the resulting .dylib libraries to the expected names (.so). For example:

.. code-block:: bash

        $ ln -s libmultinest_mpi.dylib libmultinest_mpi.so

A discussion on installing on Mac can be found in `<https://github.com/JohannesBuchner/PyMultiNest/issues/10>`_.

