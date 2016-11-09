from sys import argv, exit
from getopt import *
import numpy as np
from pbc_tools import *
from conf_tools import *
from regions import inBox, inRange
import argparse

# def cut(atoms, chainLength, box):
#   n = 0
#   molecule = []
#   drop = []
#   newBox = [[-0.5*R+comp,0.5*R+comp] for comp in origin]
#   # loop over each molecule
#   for atom in atoms:
#     if n == chainLength:
#       molecule = np.array(molecule)
#       # unnecessary because of unwrap
#       #molecule = unwrap(molecule, box)
#       # calculate its center of mass
#       COM = [calcCOM(molecule[:,d]) for d in range(len(molecule[0]))]
#       # check whether it's in the region
#       if inBox(COM, newBox): # if inRange(COM[0],[-0.5*R,0.5*R]): # inCylinder(COM, R, 0.5*R, origin):
#         # add it to drop if it is
#         drop.extend(molecule)
#       n = 0
#       molecule = []
#     n += 1
#     molecule.append(atom)
# 
#   return drop

def filter_atoms(atoms, atomtype):

  # select atoms of type atomtype
  atoms = atoms[atoms[:,2] == atomtype]

  # sort atoms
  atoms = atoms[atoms[:,0].argsort()]
  atoms[:,:2] = atoms[:,:2] + 1 - atoms[0,:2]

  return atoms

def cut(atoms, chainLength, atomtype, box):
  
  # calculate COM for each molecule
  nMol = len(atoms)/chainLength
  atoms = atoms.reshape((nMol, chainLength, 6))
  com = atoms[:,:,3:].mean(1)

  # check that each COM is in the box
  interior = ((com > box[:,0]) & (com < box[:,1]))
  interior = (interior.sum(1) == 3)

  # collect all molecules with COM inside box
  atoms = atoms[interior]

  # make the bonds
  bonds = generate_bonds(atoms, atomtype)

  # reshape the atoms array so that it has only 2 axes
  atoms = atoms.reshape((-1,6))

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
    box, atoms, bonds = readConf(inp,['1'])

  atoms = np.array(atoms,dtype=float)
  box = np.array(box,dtype=float)
  bonds = np.array(bonds,dtype=float)

  atoms = filter_atoms(atoms, atomtype)

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
             {"atoms":1,"bonds":1}, 
             [1],
             'polymer film\n' 
            )

if __name__ == "__main__":
  main()
