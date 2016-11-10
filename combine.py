from sys import argv, exit
from writeAtoms import *
from pbc_tools import *
from conf_tools import *
import numpy as np
import argparse

def getArgs(argv):

  parser = argparse.ArgumentParser(description="get filenames to combine files")
  parser.add_argument("-i", "--input", nargs=2)
  parser.add_argument("-o", "--output")
 
  return parser.parse_args()

def formatOutput(atoms, bonds, types):
  # Separates top part of trajectory
  #top = [atom[:2] + [str(int(atom[2]) + 0)] + atom[3:] for atom in atoms[:types[0]]]
  top = [atom for atom in atoms[:types[0]]]
  # Separates bottom part of trajectory
  #           str(                                                      )
  #               int(       )+int(          )-int(                  )+1
  #                   atom[0]      top[-1][0]      atoms[        ][0]
  #                                                      types[0]
  bottom  = [[str(int(atom[0])+int(top[-1][0])-int(atoms[types[0]][0])+1),
              str(int(atom[1])+int(top[-1][1])-int(atoms[types[0]][1])+1)
             ] + atom[2:] for atom in atoms[types[0]:]]
  # # fixes bonds
  # bonds = [bond[:2] + [str(int(bond[2])+int(top[-1][0])-int(atoms[types[0]][0])+1),
  #                      str(int(bond[3])+int(top[-1][0])-int(atoms[types[0]][0])+1)]
  #                      for bond in bonds] 

  return top + bottom, bonds

def makeBox(box1,box2):
  box = box1[:2]
  box.append([0.0,0.0])
  if box1[2][0] > box2[2][0]:
    box[2][0] = box2[2][0]
  else:
    box[2][0] = box1[2][0]
  if box1[2][1] < box2[2][1]:
    box[2][1] = box2[2][1]
  else:
    box[2][1] = box1[2][1]

  return box

def translateAtoms(atoms, offset):
  # select coordinates and convert to float
  coords = [[float(i) for i in atom[3:]] for atom in atoms]
  # add coordinates and offset
  coords = [[i + j for (i,j) in zip(line,offset)] for line in coords]
  # convert coordinates back to strings
  coords = [[str(i) for i in line] for line in coords]
  # recombine the lists
  return [atom[:3] + coord for (atom,coord) in zip(atoms,coords)]

def boundAtoms(atoms):
  coords = [[float(i) for i in atom[3:]] for atom in atoms]
  bounds = [coords[0][2],coords[0][2]]
  for coord in coords:
    if coord[2] < bounds[0]:
      bounds[0] = coord[2]
    elif coord[2] > bounds[1]:
      bounds[1] = coord[2]

  return bounds

def main():

  files = getArgs(argv)

  nParticles = 63

  atomSet = '3'
  with open(files.input[0], 'r') as inp1:
    box1, atoms1, bonds = readConf(inp1,atomSet)

  atomSet = '3'
  with open(files.input[1], 'r') as inp2:
    box2, atoms2, bonds = readConf(inp2,atomSet)

  # # calculates box side lengths
  # lengths = [float(ax[1])-float(ax[0]) for ax in box1]
  # # unwrap atomic coordinates (assumes that there are periodic image indices)
  # coords2 = [[float(xyz[3+i]) + float(xyz[6+i])*lengths[i] for i in range(3)] for xyz in atoms2]
  # # combines atomic information with unwrapped coordinates
  # atoms2 = [atoms2[i][:3] + [str(xyz) for xyz in coords2[i]] for i in range(len(coords2))]

  # change atom id number to type int
  atoms1 = [[int(row[0])] + row[1:6] for row in atoms1] 
  atoms2 = [[int(row[0])] + row[1:6] for row in atoms2] 

  # sort atoms
  atoms1.sort()
  atoms2.sort()

  # change atom id number back to type str
  atoms1 = [[str(row[0])] + row[1:] for row in atoms1] 
  atoms2 = [[str(row[0])] + row[1:] for row in atoms2] 

  # find the vertical bounds of the two configurations
  bounds1 = boundAtoms(atoms1)
  bounds2 = boundAtoms(atoms2)

  # shift the bottom atoms to the appropriate location and leave a small gap between confs
  margin = 0 #-1
  shift = bounds1[1]-bounds2[0]+margin
  atoms2 = translateAtoms(atoms2, [0,0,shift]) 

  # join and format the configuration information
  atoms = atoms1 + atoms2
  types = [len(atoms1),len(atoms2)]
  atoms, bonds = formatOutput(atoms, bonds, types)

  # construct the box dimensions
  box = box1[:2] + [[str(bounds2[0]+shift),str(bounds1[1])]]
  # title = 'saturated infiltration simulation with cylindrical capillary' 
  title = 'cylindrical capillary + flat base' 
  # ntypes = 3
  ntypes = 1

  write_xyz(files.output, [line[2:] for line in atoms])
  write_conf(files.output, atoms, bonds, box, {"atoms":ntypes, "bonds":0}, [1] * ntypes, title)

  # # check that the bonds and atoms correspond properly
  # test_atoms = np.array(atoms)
  # test_atoms = set(test_atoms[test_atoms[:,2] == '1'][:,0])
  # test_bonds = np.array(bonds)
  # aSet = set(test_bonds[:,2])
  # bSet = set(test_bonds[:,3])
  # bSet = bSet | aSet
  # test = test_atoms - bSet
  # print test

if __name__ == "__main__":
  main()
