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

def main():

  parser = argparse.ArgumentParser(description="build configuration by combining files")
  parser.add_argument('input')
  parser.add_argument('output')
  args = parser.parse_args()

  configurations = {'files':[],'types':[],'padding':[]}
  with open(args.input, 'r') as file:
    for line in file:
      # remove comments
      line = line.split('#')[0]
      # split entries by whitespace
      line = line.split()
      if len(line) < len(configurations)
      
        
        

  with open(args.input[0], 'r') as inp1:
    box1, atoms1, bonds1 = readConf(inp1,args.atomtypes)

  with open(args.input[1], 'r') as inp2:
    box2, atoms2, bonds2 = readConf(inp2, '1') #args.atomtypes)

  # sort atoms
  atoms1, atoms2 = atoms1[atoms1[:,0].argsort()], atoms2[atoms2[:,0].argsort()]

  # select atoms
  #atoms2 = atoms2[atoms2[:,2] == 2]
  #atoms2 = atoms2[np.logical_and(atoms2[:,3] >= box1[0,0], atoms2[:,3] < box1[0,1])]
  #atoms2 = atoms2[np.logical_and(atoms2[:,4] >= box1[1,0], atoms2[:,4] < box1[1,1])]

  # unwrap atoms if the information is present
  try:
    atoms2[:,3:6] = atoms2[:,3:6] + atoms2[:,6:]*(box2[:,1] - box2[:,0])
  except ValueError:
    pass

  # find the vertical bounds of the two configurations
  bounds1, bounds2 = boundAtoms(atoms1[:,3:]), boundAtoms(atoms2[:,3:])

  # shift the bottom atoms to the appropriate location and leave a small gap between confs
  #margin = 0 
  #shift = bounds1[1]-bounds2[0]+margin
  #atoms2[:,3:6] = atoms2[:,3:6] + np.array([0,0,shift])
  shift = 0
  atoms1[:,5] = atoms1[:,5] - bounds1[0]
  atoms2[:,5] = atoms2[:,5] - bounds2[1]

  # join and format the configuration information
  atoms = np.zeros((atoms1.shape[0]+atoms2.shape[0],9))
  atoms[:atoms1.shape[0],:6] = atoms1[:,:6]
  atoms[atoms1.shape[0]:,:6] = atoms2[:,:6]
  atoms[:,0] = np.arange(len(atoms))+1

  bonds = bonds2
  # bonds = np.zeros((bonds1.shape[0]+bonds2.shape[0],4))
  # bonds[:bonds1.shape[0]] = bonds1
  # bonds[bonds1.shape[0]:] = bonds2

  # construct the box dimensions
  box = np.zeros((3,2))
  box[:2] = box2[:2]
  box[2] = [bounds2[0]-bounds2[1],bounds1[1]-bounds1[0]+shift]
  title = 'surface and cylinder with R=11; replication of PRL/Binder paper'
  ntypes = 4

  write_xyz(args.output, [line[2:6] for line in atoms])
  write_conf(args.output, atoms, bonds, box, {"atoms":ntypes, "bonds":1}, [1] * ntypes, title)
  write_traj(args.output, np.delete(atoms[:,:6],1,1), box, mode='w')

if __name__ == "__main__":
  main()
