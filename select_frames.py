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

def write_height(output, height, mode):

  # open output file
  with open(output, mode) as file:
    # join and write height and time
    file.write(height)

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input", required=True)
  parser.add_argument("-o", "--output", required=True)
  parser.add_argument("-t", "--step", type=int)
  parser.add_argument("-n", "--frames", type=int)
  args = parser.parse_args()

  suffix = '_sparse'

  with open(args.output+suffix, 'w') as file:
    file.write('')

  with open(args.input, MODE) as file:
    for n in range(args.frames):
      # time, frame = readFrame(file)
      time, height = read_height(file)
      # time, density = read_density(file)
      if time % args.step == 0:
        # write_xyz(args.output+suffix, frame, time, 'a')
        write_height(args.output+suffix, height, 'a')
        # write_density(args.output+suffix, density, time, 'a')

if __name__ == '__main__':
  main()
