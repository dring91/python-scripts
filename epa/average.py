#!/usr/bin/python

import numpy as np

def command(data,*args):
  default = 100

  data = data.astype(float)
  try:
    data = data[data[:,0] > int(args[0])]
  except IndexError:
    data = data[data[:,0] > default]
  mean = data[:,1].mean()
  stdev = data[:,1].std()
  
  return mean, stdev
