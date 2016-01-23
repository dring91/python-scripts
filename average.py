#!/usr/bin/python

import numpy as np
from sys import argv

def main():

  for name in argv[1:]:
    data = []
    with open(name,'r') as file:
      for line in file:
        L = line.split()
        if len(L) > 0 and L[0] != '#':
          data.append(L)

    data = np.array(data,dtype=float)
    drange = (data[:,0] > 23000) & (data[:,0] < 30000)
    data = data[drange]
    print data[:,1].mean(), data[:,1].std()

if __name__ == '__main__':
  main()
