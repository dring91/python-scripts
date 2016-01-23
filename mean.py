#!/usr/bin/python

import numpy as np

def command(data, over=None):

  if over == 'rows':
    over = 0
  elif over == 'cols':
    over = 1

  data = np.array(data,dtype=float)

  return [data.mean(over)[0], np.average(data[:,0], weights=data[:,1])]
