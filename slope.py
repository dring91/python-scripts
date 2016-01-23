#!/usr/bin/python

# write a script that calculates tangent slopes
from sys import argv
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

def fitSubset(x,y):
  X = np.array([np.ones_like(x),x,x**2])
  XTX = np.einsum('ik,jk',X,X)
  XTY = np.einsum('k,jk',y,X)
  Beta = inv(XTX).dot(XTY)

  return X, Beta

def main():

  path = argv[1]
  filename = argv[2]

  # Read in height data
  data = []
  with open(path+'height_'+filename,'r') as file:
    for line in file:
      L = line.split()
      if len(L) == 0:
        break
      if len(L) > 0 and L[0] != '#':
        data.append(L)

  data = np.array(data,dtype=float)
  data[:,0] = data[:,0]*0.002/1000
  data[:,1] = data[:,1]**2
  data[:,2] = data[:,2]**2
  # data[:,3] = data[:,3]**2
  
  step = 50
  slopes = []
  for i,point in enumerate(data[:-step,:]):
    X0 = step/2+i
    (model, fit) = fitSubset(data[i:step+i,0],data[i:step+i,1])
    slope85 = fit[1]+2*fit[2]*data[X0,0]
    (model, fit) = fitSubset(data[i:step+i,0],data[i:step+i,2])
    slope99 = fit[1]+2*fit[2]*data[X0,0]
    # (model, fit) = fitSubset(data[i:step+i,0],data[i:step+i,3])
    # slopePoly = fit[1]+2*fit[2]*data[X0,0]

    slopes.append([data[X0,0],slope85,slope99]) # ,slopePoly])

  slopes = np.array(slopes)
  plt.figure()
  plt.plot(slopes[:,0],slopes[:,1],'o')
  plt.plot(slopes[:,0],slopes[:,2],'x')
  # plt.plot(slopes[:,0],slopes[:,3],'s')
  plt.figure()
  plt.plot(data[:,0],data[:,1],'o')
  plt.plot(data[:,0],data[:,2],'x')
  # plt.plot(data[:,0],data[:,3],'s')
  plt.show()

  # with open(path+'slope_'+filename,'w') as file:
  #   file.write('#  t  slope85  slope99 slopePoly\n')
  #   [file.write('%d %f %f %f\n' % tuple(line)) for line in slopes]

if __name__ == '__main__':
  main()
