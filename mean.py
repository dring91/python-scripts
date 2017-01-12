import numpy as np
import argparse

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input")
  args = parser.parse_args()

  with open(args.input, "r") as file:
    data = np.loadtxt(file)

  mean = data.mean(0)
  var = data.var(0)
  stdev = np.sqrt(var)

  with open(args.input, "a") as file:
    np.savetxt(file, mean.reshape((1,5)), "%d %f %f %f %f")
    np.savetxt(file, stdev.reshape((1,5)), "%d %f %f %f %f")

if __name__ == '__main__':
  main()
