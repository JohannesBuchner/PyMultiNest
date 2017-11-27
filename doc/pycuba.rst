PyCuba
=======================================

Cuba is a library for multidimensional numerical integration.
The `manual <https://github.com/JohannesBuchner/cuba/raw/master/cuba.pdf>`_ 
in `Cuba <http://www.feynarts.de/cuba/>`_ is helpful as a reference.
Parallelisation is currently not supported.

The algorithms available are

* Cuhre
* Divonne
* Suave
* Vegas

Make sure you `installed everything <install>`_ correctly.

PyCuba Demo Program
----------------------
Take a look at the demo

- `Demo pycuba_demo.py <https://github.com/JohannesBuchner/PyMultiNest/blob/master/pycuba_demo.py>`_.


PyCuba API
------------------

The return value is always a dictionary:

.. code-block:: python

    { 'neval': number of evaluations, 
      'fail': 0 or 1,
      'comp': number of components, usually 1,
      'nregions': number of regions used,
      'results': [ # a list of results for each component
         {
           'integral': value of the integral, 
           'error':  uncertainty, 
           'prob': probability (see manual),
         }, 
         ... 
      ]
    }

Define your integrand, like so:

.. code-block:: python

	def Integrand(ndim, xx, ncomp, ff, userdata):
		# access the current parameters
		x,y,z = [xx[i] for i in range(ndim.contents.value)]
		# compute the result
		result = math.sin(x)*math.cos(y)*math.exp(z)
		# store the result (here only one component)
		ff[0] = result
		return 0

It will be called with xx in the interval from 0 to 1 (so scale to your borders
in this function).

The :py:func:`pycuba.demo` in the pycuba source shows how to run all of the algorithms and access the results:

.. literalinclude:: ../pycuba/__init__.py
	:language: python
	:start-after: def demo():
	:end-before: if __name__ == '__main__':

API documentation
------------------------

.. automodule:: pycuba
	:members:

