from sys import argv
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
import argparse

def fitSubset(x,y):
  # model to fit
  X = np.array([np.ones_like(x),x,x**2])
  #X = np.array([np.ones_like(x),x])
  # linear solution to model
  XTX = np.einsum('ik,jk',X,X)
  XTY = np.einsum('k,jk',y,X)
  Beta = inv(XTX).dot(XTY)

  return X, Beta

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('-p','--path')
  parser.add_argument('-i','--input')
  parser.add_argument('-o','--output')
  args = parser.parse_args()

  # Read in height data
  data = []
  with open(args.path+'height_'+args.input,'r') as file:
    for line in file:
      L = line.split()
      if len(L) == 0:
        break
      if len(L) > 0 and L[0] != '#':
        data.append(L)

  data = np.array(data,dtype=float)
  #data[:,0] = data[:,0]*0.002/1000
  data[:,0] = np.log10(data[:,0])
  data[:,3] = np.log10(data[:,3])
  data = data[1:]
  
  step = 10
  slopes = []
  for i,point in enumerate(data[:-step,:]):
    X0 = step/2+i
    (model, fit) = fitSubset(data[i:step+i,0],data[i:step+i,3])
    prediction = np.einsum('i...,i',model,fit)
    slope = fit[1]+2*fit[2]*data[X0,0]
    #slope = fit[1]*data[X0,0]

    #plt.figure()
    #plt.plot(data[i:step+i,0],prediction,'o')
    #plt.plot(data[i:step+i,0],data[i:step+i,3],'o')
    #plt.show()

    slopes.append([data[X0,0],slope])

  slopes = np.array(slopes)
  plt.figure()
  plt.plot(slopes[:,0],slopes[:,1],'o')
  #plt.figure()
  plt.plot(data[:,0],data[:,3],'o')
  plt.show()

  # with open(path+'slope_'+filename,'w') as file:
  #   file.write('#  t  slope85  slope99 slopePoly\n')
  #   [file.write('%d %f %f %f\n' % tuple(line)) for line in slopes]

if __name__ == '__main__':
  main()
