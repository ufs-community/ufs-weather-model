#!/usr/bin/env python
import sys 
import numpy as np
from netCDF4 import Dataset

with Dataset(sys.argv[1]) as nc1, Dataset(sys.argv[2]) as nc2:
  for varname in nc1.variables.keys():
    diff = nc2[varname][:]-nc1[varname][:]
    if (np.abs(diff)).max() != 0:
      print("The two nc files are different.")
      sys.exit(1)
