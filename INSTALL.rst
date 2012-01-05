
1. build MultiNest in some directory $MULTINEST
  
  $ export MULTINEST=/my/multinest/directory

  $ cd $MULTINEST && make libnest3.so

1A. build APEMoST in some directory $APEMOST
  
  $ export APEMOST=/my/apemost/directory

  $ cd $APEMOST && make libapemost.so

2. go into multinest_bridge and build 

  $ cd multinest_bridge && make libcnest.so

3. put all of these directories in the dynamic library load path

  $ export LD_LIBRARY_PATH=$MULTINEST:$PWD/multinest_bridge:$APEMOST

4. tell python where to find the libraries

  $ export PYTHONPATH=$PWD

5. go to some directory and run your stuff

  $ cd /tmp
  $ python -c 'import pymultinest'
  $ python -c 'import pyapemost'
  $ python -c 'import pymultinest_demo'
  $ python -c 'import pyapemost_demo'


