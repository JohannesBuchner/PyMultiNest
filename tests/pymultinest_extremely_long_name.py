# Script to test MultiNest failure when running with an extremely long file name (> than 1000 characters)

import pytest
from pymultinest.solve import solve

def test():
    with pytest.raises(ValueError):
        solve(lambda cube: -0.5 *(cube**2).sum(), lambda cube: cube, 2, 
              resume=True, verbose=True,
              outputfiles_basename=1100*"a")

