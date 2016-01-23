#!/usr/bin/python

# write a script that calculates tangent slopes
from sys import argv
import numpy as np
# from numpy.linalg import inv
import matplotlib.pyplot as plt

def readFile(filename, nFrames):
  nSlices = 200
  label = np.zeros((nSlices, nFrames))
  data = np.zeros((nSlices, nFrames))
  with open(filename,'r') as file:
    n = 0
    for line in file:
      L = line.split()
      if len(L) > 1 and L[0] != '#':
        i = n % nSlices
        j = n / nSlices
        label[i,j] = float(L[0])
        data[i,j] = float(L[1])
        n += 1
  # data.astype(float)
  # data = np.array(data, dtype=float)

  return label, data

def main():
  # load data
  porosityFile = argv[1]
  densityFile = argv[2]

  zi, porosity = readFile(porosityFile, 12)
  zj, density = readFile(densityFile,366)
  # print density.shape

  # filter density data
    # bulkDensity?
    # 200 points with greatest z value
  density = density[:,range(0,360,30)]
  zj = zj[:,range(0,360,30)]
  for li,lj,coli,colj in zip(zi,zj,porosity, density):
    colj = colj[lj > 0]
    for j in colj:
      pass
      # check for j <= i

  # # build integrand
  # integrand = density / porosity
  # # choose integration method
  # # perform integration
  #   # for this case, a simple average should do
  #   average = integrand.mean() / bulkDensity
  # # output results
  # print average

if __name__ == '__main__':
  main()
