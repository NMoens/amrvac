#!/usr/bin/env python

# Purpose: return 2-norm, 1-norm and inf-norm between two VTK files
#
# This script was inspired by script at tests/rho/convergence under the 'old'
# amrvac branch
#
# Author: Jannis Teunissen

import argparse
import vtk as v
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import sys


def get_args():
    # Get and parse the command line arguments
    p = argparse.ArgumentParser(
        description='Return 2-norm, 1-norm and inf-norm between two files')
    p.add_argument('file_a', type=str, help='First VTK file')
    p.add_argument('file_b', type=str, help='Second VTK file')
    p.add_argument('var', type=str, help='Variable name')
    p.add_argument('-vtktype', type=str, default='unstructured',
                   choices=['unstructured'], help='Type of VTK data')
    return p.parse_args()


args = get_args()

if args.vtktype == 'unstructured':
    datareader = v.vtkXMLUnstructuredGridReader()
else:
    print('Cannot handle this type of VTK file')
    sys.exit(1)

datareader.SetFileName(args.file_a)
datareader.Update()
d_a = datareader.GetOutput()
var_a = vtk_to_numpy(d_a.GetCellData().GetArray(args.var))

datareader.SetFileName(args.file_b)
datareader.Update()
d_b = datareader.GetOutput()
var_b = vtk_to_numpy(d_b.GetCellData().GetArray(args.var))

# Output 2-norm 1-norm inf-norm
print(np.linalg.norm(var_a - var_b),
      np.linalg.norm(var_a - var_b, 1),
      np.linalg.norm(var_a - var_b, np.inf))
