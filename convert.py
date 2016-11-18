import numpy as np
import argparse
from conf_tools import *

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input")
  parser.add_argument("-o","--output")
  parser.add_argument("-n","--frames", type=int)
  args = parser.parse_args()

  with open(args.output+'.lammpstrj', "w") as file:
    file.write('')

  with open(args.input, "r") as file:
    for n in range(args.frames):
      time, atoms = readFrame(file)

      box = np.zeros((3,2))
      box[:,0] = atoms[:,1:].min(axis=0)
      box[:,1] = atoms[:,1:].max(axis=0)

      data = np.zeros((atoms.shape[0],5))
      data[:,0] = np.arange(atoms.shape[0])+1
      data[:,1] = atoms[:,0]
      data[:,2:] = (atoms[:,1:] - box[:,0])/(box[:,1] - box[:,0])

      write_traj(args.output, data, box, time)

if __name__ == '__main__':
  main()
