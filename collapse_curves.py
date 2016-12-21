from sys import argv,exit
from conf_tools import *
from pbc_tools import *
import argparse

MODE = 'r'

def read_density(file):

  # define the dimensions of the density histogram matrix
  nBins = 152
  nCols = 4

  # read header and extract time
  L = file.readline().split()
  time = int(L[-1])
  file.readline()

  # read density histogram
  density = np.fromfile(file, float, nBins*nCols, " ")
  density = density.reshape((nBins, nCols))

  return time, density

def write_density(output, density, time, mode):

  with open(output, mode) as file:
    file.write('# %d\n' % time)
    np.savetxt(file, density, "%.8f")
    file.write('\n\n')

def read_height(file):

  # read line in from height file
  line = file.readline()
  # check that the line is non-empty and not a comment
  if len(line) > 0: 
    while line[0] == "#":
      line = file.readline()
    L = line.split()
    time = int(L[0])

  return time, line

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument("-z", "--height", required=True)
  parser.add_argument("-d", "--density", required=True)
  parser.add_argument("-o", "--output", required=True)
  parser.add_argument("-n", "--frames", type=int)
  args = parser.parse_args()

  with open(args.output, 'w') as file:
    file.write('')

  with open(args.height, MODE) as hfile, open(args.density, MODE) as dfile:
    for n in range(args.frames):
      time_h, height = read_height(hfile)
      time_d, density = read_density(dfile)
      
      if time_h == time_d:
      #write_density(args.output+suffix, density, time, 'a')

if __name__ == '__main__':
  main()
