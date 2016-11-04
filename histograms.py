import numpy as np
from numpy.random import rand, randn
import matplotlib.pyplot as plt

def main():
  
  # generate random data to test method
  indices = np.arange(200)
  data = rand(200,2)
  # transform data using a function to test how well the method captures different gradients
  line = lambda x: 0.25*x+1
  # plt.plot(data[:,0], data[:,1], 'o')
  plt.plot(randn(200)+line(indices),'o')
  plt.show()
  
  # divide data into N bins with p points inside each bin
  
  # calculate the appropriate variance of each bin
  # split bins with 'high' variance and recalculate variances in modified bins

if __name__ == '__main__':
  main()
