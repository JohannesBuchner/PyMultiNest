Installing PyMultiNest and PyAPEMoST
=================================================

Installing the python library
------------------------------------

Just run ::
   
     $ python setup.py install

or if you only want to install locally::

     $ python setup.py install --user

Installing the bridge to MultiNest
------------------------------------

1. build MultiNest in some directory $MULTINEST::
  
     $ export MULTINEST=/my/multinest/directory

     $ cd $MULTINEST && make libnest3.so

2. go into multinest_bridge and build::

     $ cd multinest_bridge && make libcnest.so

Installing the bridge to APEMoST
------------------------------------

3. build APEMoST in some directory $APEMOST::
  
     $ export APEMOST=/my/apemost/directory

     $ cd $APEMOST && make libapemost.so

Running some code
--------------------------

PyMultiNest / PyAPEMoST have to be able to find the corresponding 
libraries. So put the three directories in the dynamic library load path::

     $ export LD_LIBRARY_PATH=$MULTINEST:$PWD/multinest_bridge:$APEMOST

Test importing the libraries::

     $ python -c 'import pymultinest'
     $ python -c 'import pyapemost'

5. go to some directory and run your stuff::

     $ cd /tmp
     $ python $OLDPWD/pymultinest_demo.py
     $ python $OLDPWD/pyapemost_demo.py

Generating the documentation
----------------------------

Go into doc and run make::

     $ cd doc && make html


