from sys import argv, exit
import numpy as np
import argparse

from pbc_tools import *
from conf_tools import *

def boundAtoms(atoms):
  bounds = np.zeros(2)
  bounds[0] = atoms[:,2].min()
  bounds[1] = atoms[:,2].max()

  return bounds

def getArgs(argv):

  parser = argparse.ArgumentParser(description="combine configuration files")
  parser.add_argument("-i", "--input", nargs=2)
  parser.add_argument("-o", "--output")
  parser.add_argument("-a", "--atomtypes", nargs='+')
 
  return parser.parse_args()

def main():

  args = getArgs(argv)

  with open(args.input[0], 'r') as inp1:
    box1, atoms1, bonds1 = readConf(inp1,args.atomtypes)

  with open(args.input[1], 'r') as inp2:
    box2, atoms2, bonds2 = readConf(inp2,args.atomtypes)

  # sort atoms
  atoms1 = atoms1[atoms1[:,0].argsort()]
  atoms2 = atoms2[atoms2[:,0].argsort()]

  # unwrap atoms if the information is present
  atoms1[:,3:6] = atoms1[:,3:6] + atoms1[:,6:]*(box1[:,1] - box1[:,0])

  # find the vertical bounds of the two configurations
  bounds1 = boundAtoms(atoms1[:,3:])
  bounds2 = boundAtoms(atoms2[:,3:])

  # shift the bottom atoms to the appropriate location and leave a small gap between confs
  margin = 0 
  shift = bounds1[1]-bounds2[0]+margin
  atoms2[:,3:6] = atoms2[:,3:6] + np.array([0,0,shift])

  # join and format the configuration information
  atoms = np.zeros((atoms1.shape[0]+atoms2.shape[0],6))
  atoms[:atoms1.shape[0]] = atoms1[:,:6]
  atoms[atoms1.shape[0]:] = atoms2[:,:6]
  bonds = bonds1
  # bonds = np.zeros((bonds1.shape[0]+bonds2.shape[0],4))
  # bonds[:bonds1.shape[0]] = bonds1
  # bonds[bonds1.shape[0]:] = bonds2

  # construct the box dimensions
  box = np.zeros((3,2))
  box[:2] = box1[:2]
  box[2] = [bounds1[0],bounds2[1]+shift]
  title = 'cylindrical capillary + flat base' 
  ntypes = 3

  write_xyz(args.output, [line[2:] for line in atoms])
  write_conf(args.output, atoms, bonds, box, {"atoms":ntypes, "bonds":1}, [1] * ntypes, title)
  # write_traj(args.output, np.delete(atoms,1,1), box)

if __name__ == "__main__":
  main()
