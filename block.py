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

  sample = sample[:,1]
  sample = sample.reshape((200,5))
  mean = sample.mean(1)
  std = sample.std(1)
  # print mean, std

  plt.errorbar(actual[:,0],mean,yerr=std,fmt='o')
  plt.plot(actual[:,0],actual[:,1])
  plt.show()

if __name__ == '__main__':
  main()
