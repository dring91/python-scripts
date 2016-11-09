from sys import argv, exit
from getopt import *
import numpy as np
from pbc_tools import *
from conf_tools import *
from regions import inBox, inRange
import argparse

def cut(atoms, chainLength, box):
  n = 0
  molecule = []
  drop = []
  newBox = [[-0.5*R+comp,0.5*R+comp] for comp in origin]
  # loop over each molecule
  for atom in atoms:
    if n == chainLength:
      # unnecessary because of unwrap
      molecule = np.array(molecule)
      molecule = unwrap(molecule, box)
      # calculate its center of mass
      COM = [calcCOM(molecule[:,d]) for d in range(len(molecule[0]))]
      # check whether it's in the region
      if inBox(COM, newBox): # if inRange(COM[0],[-0.5*R,0.5*R]): # inCylinder(COM, R, 0.5*R, origin):
        # add it to drop if it is
        drop.extend(molecule)
      n = 0
      molecule = []
    n += 1
    molecule.append(atom)

  return drop

def main():

  # parse command line arguments
  parser = argparse.ArgumentParser(description='read information for script exec')
  parser.add_argument('-i', '--input', required=True)
  parser.add_argument('-o', '--output', required=True)
  parser.add_argument('-l', '--length', type=int, required=True)
  parser.add_argument('-x', nargs='2', type=int)
  parser.add_argument('-y', nargs='2', type=int)
  parser.add_argument('-z', nargs='2', type=int)
  args = parser.parse_args()

  # exit if no x, y, or z dimensions provided
  if len(args) == 3:
    print "no change to initial configuration"
    print "no file written"
    sys.exit()

  # read in configuration file
  with open(args.input, 'r') as inp:
    box, atoms, bonds = readConf(inp,['1'])

  atoms = np.array(atoms,dtype=float)
  box = np.array(box,dtype=float)
  bonds = np.array(bonds,dtype=float)

  # reconstruct box using the new dimensions provided
  try:
    box[0,:] = args.x
  except AttributeError:
    pass

  drop = cut(atoms, args.length, box)

  # generate a new set of bonds based on the new set of atoms
  bonds = np.array([])

  # write_xyz(filename, atoms, {time, mode})
  write_xyz(args.output, drop)
  # write_conf(filename, atoms, {bonds,box,types,masses,title}):
  write_conf(args.output, 
             drop,
             bonds,
             box,
             {"atoms":1,"bonds":1}, 
             [1],
             'polymer film with dimensions \n' 
            )

if __name__ == "__main__":
  main()
