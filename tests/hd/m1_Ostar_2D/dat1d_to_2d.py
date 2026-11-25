#!/usr/bin/env python3
"""
1D -> 2D spherical snapshot conversion; kept for existing commands.
This is dat_lift.py with the target dimensionality fixed to 2 (same options,
see `python3 dat_lift.py -h`). For 1D -> 3D or 2D -> 3D use dat_lift.py --ndim 3.

Example
  python3 dat1d_to_2d.py ../m1_Ostar_1D/resscan/n1280_0100.dat --reset \\
          --par parfiles_hpc/parfiles/base.par parfiles_hpc/parfiles/Star8.par
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dat_lift import main  # noqa: E402

if __name__ == '__main__':
    main(default_ndim=2)
