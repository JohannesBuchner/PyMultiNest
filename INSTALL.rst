Installing PyMultiNest and PyAPEMoST
=================================================

.. contents::

To use PyMultiNest or PyAPEMoST, follow these installation steps.
If you only want to use one of them, skip steps according.

Install on Linux
-----------------------------------

  * If you want to use PyMultiNest: 
    Get and compile MultiNest (use the cmake version from `here <https://github.com/JohannesBuchner/MultiNest>`_). You will have a lib/libmultinest.so. 
    Include the containing directory in your LD_LIBRARY_PATH.
  
  * If you want to use PyAPEMoST: 
    Get and compile `APEMoST <http://apemost.sourceforge.net/doc/>`_. You will have a libapemost.so file.
    Include the containing directory your LD_LIBRARY_PATH.
  
  * If you want to use PyCuba:
    Get and compile `Cuba <http://www.feynarts.de/cuba/>`_. You will have a libcuba.so file.
    Include the containing directory your LD_LIBRARY_PATH.
    
  * Install the python library from PyPI::

  	$ pip install pymultinest
    
    Alternatively, download the PyMultiNest source and install it locally by
    running the following commands in the PyMultiNest directory::
    
    	$ pip install ctypesGsl # a prerequisite
    	$ python setup.py install
    
    If you only want to install locally, use::
    
    	$ pip install --user ctypesGsl # a prerequisite
    	$ python setup.py install --user

Running some code
--------------------------

PyMultiNest, PyAPEMoST and PyCuba have to be able to find the corresponding 
libraries. So put the three directories in the dynamic library load path::

     $ export LD_LIBRARY_PATH=/my/multinest/directory/lib:/my/cuba/directory/:/my/apemost/directory/:$LD_LIBRARY_PATH

Test importing the libraries::

     $ python -c 'import pymultinest'
     $ python -c 'import pyapemost'

Try out the demo programs distributed in the package. They produce a lot of output, so lets create a temporary directory::

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

Go into doc and run make::

     $ cd doc && make html

Point your web browser to _build/html/index.html in doc.


Install on Mac
-----------------------------
Follow the above instructions. Use DYLD_LIBRARY_PATH instead of LD_LIBRARY_PATH.
rename the resulting .dylib libraries to the expected names (.so). For example ::

        $ ln -s libmultinest_mpi.dylib libmultinest_mpi.so

A discussion on installing on Mac can be found in `<https://github.com/JohannesBuchner/PyMultiNest/issues/10>`_.

