from sys import argv, exit
import numpy as np
import argparse

from pbc_tools import *
from conf_tools import *

def main():

  parser = argparse.ArgumentParser(description="combine configuration files")
  parser.add_argument("-i", "--input", nargs='+',
                      help="input files to combine")
  parser.add_argument("-l", "--layers", nargs='+', type=int,
                      help="the order in which to combine layers; layers are combined \
                            according to the position of their input file")
  parser.add_argument("-z", "--zero", nargs=2, type=int,
                      help="the layer and position of zero \
                            (eg - 1 0 is the bottom (0) of the second layer (1)")
  parser.add_argument("--pad", action='store_true', help="pad between layers")
  parser.add_argument("-o", "--output", help="output file name")
  parser.add_argument("-d", "--descrip", help="description for output file header")
  parser.add_argument("-f", "--formats", nargs='+', help="output formats",
                      choices=['xyz','lammps','conf'], default='conf')
  parser.add_argument("--no-image-flags",dest="entries",action="store_const",const=6,
                      default=9,help="switch to adjust number of entries for data flags")
  parser.add_argument("--use-largest-box",dest="size",action="store_const",const=np.argmax,
                      default=np.argmin,help="switch to specify how to choose box sizes")
  args = parser.parse_args()

  # initialize list of layers
  layers = []
  bounds = []
  boxes = []
  bonded = []
  # read in all input files
  for name in args.input:
    with open(name, 'r') as file:
      box, atoms, bonds = readConf(file)
      # sort atoms and bonds
      atoms = atoms[atoms[:,0].argsort()]
      if bonds.shape[0] != 0: bonds = bonds[bonds[:,2].argsort()]
      # create a set of all atomtypes
      types, mask = np.unique(atoms[:,2],return_index=True)
      # loop through each atomtype and separate out each layer
      for type in atoms[mask][:,2]:
        layer = atoms[atoms[:,2] == type] #[:,:6]
        layers.append(layer)
        bound = [min(layer[:,5]), max(layer[:,5])]
        bounds.append(bound)
        if len(bonds) > 0:
          bonded.append(bonds[bonds[:,1] == type])
        else:
          bonded.append(np.array([]))
      boxes.append(box[:2])

  # select layers and reorder
  layers = np.array(layers)
  bounds = np.array(bounds)
  layers = layers[args.layers]
  bounds = bounds[args.layers]
  bondtypes = []
  for i,(layer,bonds) in enumerate(zip(layers,bonded)):
    layer[:,2] = i+1
    if bonds.shape[0] != 0: 
      bonds[:,1] = i+1
      bondtypes.append(i+1)
  thicknesses = bounds[:,1] - bounds[:,0]

  # pad layers
  padding = int(args.pad) * np.ones_like(args.layers)
  padding[1] = 1
  thicknesses += padding

  # Calculate appropriate box size
  boxes = np.array(boxes)
  box = np.zeros((3,2))
  box[:2] = boxes[args.size(np.abs(boxes), axis=0)][0,0] #args.size(boxes,key=np.any)
  box[2] = [0,sum(thicknesses)]

  # extract zero
  l,p = tuple(args.zero)

  # apply shifts to atoms
  positions = np.cumsum(thicknesses)
  positions = np.roll(positions,1)
  positions[0] = 0
  for i,bound in enumerate(bounds):
    layers[i][:,5] += positions[i]
    layers[i][:,5] -= bound[0]
    layers[i][:,5] -= positions[l] + thicknesses[l] * p
    layers[i] = layers[i][:,:args.entries]

  # prepare atoms and bonds for writing to output
  atoms = np.concatenate(layers, axis=0)
  #bonds = np.concatenate(bonded, axis=0)
  for bond in bonded:
    if len(bond) > 0:
      bonds = np.copy(bond)
  atoms[:,0] = np.arange(atoms.shape[0])+1

  # write output files
  atomtypes = len(layers)
  bondtypes = len(bondtypes)
  if 'xyz' in args.formats:
    write_xyz(args.output, [line[2:6] for line in atoms])
  if 'conf' in args.formats:
    write_conf(args.output, atoms, bonds, box, {"atoms":atomtypes, "bonds":bondtypes}, [1] * atomtypes, args.descrip)
  if 'lammps' in args.formats:
    write_traj(args.output, np.delete(atoms[:,:6],1,1), box, mode='w')

if __name__ == "__main__":
  main()
