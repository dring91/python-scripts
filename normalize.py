#!/usr/bin/python

import numpy as np

def command(data, over=None):

  if over == 'rows':
    over = 0
  elif over == 'cols':
    over = 1

  data = np.array(data,dtype=float)

  summation = data.sum(over)*(data[1,0]-data[0,0])

  data[:,1] /= summation[1]

  return data
