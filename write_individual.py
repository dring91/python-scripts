import numpy as np
import argparse
from conf_tools import readTrj

def main():
  """ script to write each frame of a trajectory to a different file """

  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input")
  args = parser.parse_args()

  # output file information
  output = args.input.split('.')

  with open(args.input, "r") as file:
    for s in (args.steps):
      time, box, atoms = readTrj(file)

      # subsitute file extension with the time
      output[-1] = str(time)

      # sort atoms according to ID
      atoms = atoms[np.argsort(atoms[:,0])]
      # extract polymers and walls
      polymers = atoms[atoms[:,1] == 1][:,2:]
      cylinder = atoms[atoms[:,1] >= 3][:,2:]
      # change vertical origin
      polymers[:,2] -= cylinder[:,2].min()
      # select only polymers above the origin
      polymers = polymers[polymers[:,2] >= 0]

      ## write to file
      filename = '.'.join(output)
      np.savetxt(filename,polymers)

if __name__ == '__main__':
  main()
