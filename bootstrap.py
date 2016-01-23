#!/usr/bin/python

# write a script that calculates tangent slopes
from sys import argv
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

def error(actual, data):
  error = np.sqrt(np.sum((actual-data)**2)/len(data))

  return error

def main():
  sampleName = argv[1]
  actualName = argv[2]

  sample = []
  actual = []
  nPoints = 200
  nResamples = 50

  with open(sampleName,'r') as file:
    for line in file:
      L = line.split()
      if len(L) > 0 and L[0] != '#':
        sample.append(L)
  sample = np.array(sample, dtype=float)

  with open(actualName,'r') as file:
    for line in file:
      L = line.split()
      if len(L) > 0 and L[0] != '#':
        actual.append(L)
  actual = np.array(actual, dtype=float)

  # print error(actual, sample)

  # define a re-sample size
  size = 100
  indices = np.arange(size)

  bootstrap = []
  for i in range(nResamples):
    # perform random resampling on the data
    mask = np.random.choice(indices, size, replace=True)
    resample = sample[mask]
    # calculate a statistic based on the resample
    MSE = error(actual[mask][:,1], resample[:,1])
    # Add the statistic to a list that will be turned into a distribution
    bootstrap.append(MSE)

  plt.plot(bootstrap, 'o')
  plt.show()

if __name__ == '__main__':
  main()
