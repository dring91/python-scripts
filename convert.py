import numpy as np
from sys import exit
import argparse
from conf_tools import *

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input")
  parser.add_argument("-o","--output")
  parser.add_argument("-n","--frames", type=int, 
                      help="The number of frames to use; alternatively, \
                            the actual frames to use")
  args = parser.parse_args()

  ## use the file extensions to figure out the formats to use
  in_format = args.input.split('.')[-1]
  out_format = args.output.split('.')[-1]

  with open(args.output, "w") as file:
    file.write('')

  with open(args.input, "r") as file:
    for n in range(args.frames):
      if in_format == 'conf':
        time, box, atoms = readFrame(file)
      else:
        print 'input file type not supported'
        exit(1)

      box = np.zeros((3,2))
      box[:,0] = atoms[:,1:].min(axis=0)
      box[:,1] = atoms[:,1:].max(axis=0)

      data = np.zeros((atoms.shape[0],5))
      data[:,0] = np.arange(atoms.shape[0])+1
      data[:,1] = atoms[:,0]
      data[:,2:] = (atoms[:,1:] - box[:,0])/(box[:,1] - box[:,0])

      if out_format == 'lammps':
        write_traj(args.output, data, box, time)
      else:
        print 'output file type not supported'
        exit(2)

if __name__ == '__main__':
  main()
