# Script to test MultiNest when running with a long file name (> than 100 characters)

if __name__ == '__main__':
    # avoid running in pytest
    from pymultinest.solve import solve

    solve(lambda cube: -0.5 *(cube**2).sum(), lambda cube: cube, 2, 
          resume=True, verbose=True,
          outputfiles_basename="a_really_really_really_really_really_really_really_really_really_really_really_really_really_really_really_really_really_long_name-")
