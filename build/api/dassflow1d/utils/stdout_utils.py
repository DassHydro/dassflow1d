# coding: utf-8
from __future__ import print_function

import sys


def stdout_progress(prefix, value, minval=0, maxval=100, nsteps=80, symbol="#"):
  
  if value == minval:
    
    sys.stdout.write(prefix)
    sys.stdout.flush()
    
  elif int((value-1) * nsteps / (maxval - minval)) < int(value * nsteps / (maxval - minval)):
    
    nsub = int(value * nsteps / (maxval - minval)) - int((value-1) * nsteps / (maxval - minval))
    for isub in range(0, nsub):
      sys.stdout.write("%s" % symbol)
    sys.stdout.flush()
    
  if value >= maxval:
    
    sys.stdout.write("\n")
    sys.stdout.flush()
