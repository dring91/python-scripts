from sys import argv
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
import argparse

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
  data[:,0] = data[:,0]*0.002/1000
  #data[:,0] = np.log10(data[:,0])
  #data[:,3] = np.log10(data[:,3])
  #data = data[1:]

  step = (data[1:,0] - data[:-1,0]).mean()

  slopes = np.gradient(data[:,3], step)
  
  plt.figure()
  plt.plot(data[:,0],slopes,'o')
  #plt.figure()
  #plt.plot(data[:,0],data[:,3],'o')
  plt.show()

  # with open(path+'slope_'+filename,'w') as file:
  #   file.write('#  t  slope85  slope99 slopePoly\n')
  #   [file.write('%d %f %f %f\n' % tuple(line)) for line in slopes]

if __name__ == '__main__':
  main()
