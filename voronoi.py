import numpy as np
import argparse
import matplotlib.pyplot as plt
import subprocess as sproc

from conf_tools import readTrj

def main():

  """ a script to calculate the density with voronoi diagrams """

  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input")
  parser.add_argument("-s", "--steps")
  parser.add_argument("-r", "--radius", type=float, 
                      help='radius of wall used in voronoi computation')
  parser.add_argument("--from-center", action='store_true', dest='correction',
                      help='option to correct radius for particle size')
  parser.add_argument("voro", nargs=argparse.REMAINDER, 
                      help="arguments to be dispatched to voro++ commandline utility")
  args = parser.parse_args()

  with open(args.input, "r") as file:
    for s in (args.steps):
      time, box, atoms = readTrj(file)

      # sort atoms according to ID
      atoms = atoms[np.argsort(atoms[:,0])]
      # extract polymers and walls
      polymers = atoms[atoms[:,1] == 1][:,2:]
      cylinder = atoms[atoms[:,1] >= 3][:,2:]
      # change vertical origin
      polymers[:,2] -= cylinder[:,2].min()
      # select only polymers above the origin
      polymers = polymers[polymers[:,2] >= 0]

      # call voro++ commandline utility with subprocess module
      # setup shell script
      args = ['time','voro++'] + args.voro
      sproc.call(args)


if __name__ == '__main__':
  main()
