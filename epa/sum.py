#!/usr/bin/python

import numpy as np

def command(data,*args):
  # how do you generalize sum?
  # how do you incorporate weights?
  #  - ui = user weights, wi = data weights
  wi = data[:,-1]
  for arg in args:
    if arg == 'equal':
      wi = np.ones(len(data))
    if arg == 'average':
      wi = np.ones(len(data), dtype=float) / len(data)
  sum = np.dot(data.T, wi)
  
  return sum
