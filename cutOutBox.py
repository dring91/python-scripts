from sys import argv, exit
from getopt import *
import numpy as np
from pbc_tools import *
from conf_tools import *
from regions import inBox, inRange
import argparse

def filter_atoms(atoms, atomtype):

  # sort atoms
  atoms = atoms[atoms[:,0].argsort()]
  atoms[:,:2] = atoms[:,:2] + 1 - atoms[0,:2]

  # select atoms of type atomtype
  polymer = atoms[atoms[:,2] == atomtype]
  substrate = atoms[atoms[:,2] != atomtype]

  return polymer, substrate

def cut(atoms, chainLength, atomtype, box):
  
  # separate atoms by type
  polymer, substrate = filter_atoms(atoms, atomtype)
  
  # filter substrate atoms
  interior = ((substrate[:,3:5] > box[:2,0]) & (substrate[:,3:5] < box[:2,1]))
  interior = (interior.sum(1) == 2)
  substrate = substrate[interior]

  # calculate COM for each polymer molecule
  nMol = len(polymer)/chainLength
  polymer = polymer.reshape((nMol, chainLength, 6))
  com = polymer[:,:,3:].mean(1)

  # check that each COM is in the box
  interior = ((com > box[:,0]) & (com < box[:,1]))
  interior = (interior.sum(1) == 3)

  # collect all polymer molecules with COM inside box
  polymer = polymer[interior]

  # make the bonds
  bonds = generate_bonds(polymer, atomtype)

  # reshape the atoms array so that it has only 2 axes
  polymer = polymer.reshape((-1,6))

  # join the two atom arrays along their outer axis
  atoms = np.concatenate((polymer,substrate), axis=0)

  return atoms, bonds

def generate_bonds(atoms, atomtype):
  
  # generate the bonds set
  bonds = np.zeros((atoms.shape[0], atoms.shape[1]-1, 4))
  bonds[:,:,1] = atomtype
  bonds[:,:,2] = atoms[:,:-1,0]
  bonds[:,:,3] = atoms[:,1:,0]

  # reshape bonds
  bonds = bonds.reshape((-1,4))

  # number bonds
  bonds[:,0] = np.arange(atoms.shape[0]*(atoms.shape[1]-1))+1

  return bonds

def main():

  # parse command line arguments
  parser = argparse.ArgumentParser(description='read information for script exec')
  parser.add_argument('-i', '--input', required=True)
  parser.add_argument('-o', '--output', required=True)
  parser.add_argument('-l', '--length', type=int, required=True)
  parser.add_argument('-x', nargs=2, metavar=('xlo','xhi'), type=int)
  parser.add_argument('-y', nargs=2, metavar=('ylo','yhi'), type=int)
  parser.add_argument('-z', nargs=2, metavar=('zlo','zhi'), type=int)
  args = parser.parse_args()

  if (args.x, args.y, args.z) is (None, None, None):
    print "Box dimensions unchanged"
    print "New film not generated"
    sys.exit()

  atomtype = 1

  # read in configuration file
  with open(args.input, 'r') as inp:
    box, atoms, bonds = readConf(inp,['1','2'])

  atoms = np.array(atoms,dtype=float)
  box = np.array(box,dtype=float)
  bonds = np.array(bonds,dtype=float)

  # set the vertical zero
  atoms[:,5] -= box[2,0]
  print box[2,0]

  # reconstruct box using the new dimensions provided
  if args.x is not None:
    box[0,:] = args.x
  if args.y is not None:
    box[1,:] = args.y
  if args.z is not None:
    box[2,:] = args.z

  drop, bonds = cut(atoms, args.length, atomtype, box)

  # write_xyz(filename, atoms, {time, mode})
  write_xyz(args.output, drop[:,2:])
  # write_conf(filename, atoms, {bonds,box,types,masses,title}):
  write_conf(args.output, 
             drop,
             bonds,
             box,
             {"atoms":2,"bonds":1}, 
             [1,1],
             'polymer film' 
            )

if __name__ == "__main__":
  main()
