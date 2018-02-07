import numpy as np
from argparse import ArgumentParser
from conf_tools import *

def between(values, bounds):
  minima, maxima = values < bounds[0], values > bounds[1]
  return values * np.logical_or(minima,maxima)

def main():
  parser = ArgumentParser()
  parser.add_argument('-i','--input')
  parser.add_argument('-o','--output')
  parser.add_argument('-t','--types',default=[],nargs='+',type=int)
  parser.add_argument('-d','--description',default='')
  parser.add_argument('--region',action='append',nargs='+')
  args = parser.parse_args()

  with open(args.input) as file:
    types, box, atoms, bonds = readConf(file)
  # if no types are given, all types are removed from the region
  if args.types == []: args.types = types['atoms']
  # label dimensions
  dims = {'x':0,'y':1,'z':2}

  #create subparser to handle region data
  subpar = ArgumentParser()
  subpar.add_argument('dim',choices={'x','y','z'})
  subpar.add_argument('bounds',nargs=2,type=float)
  for region in map(subpar.parse_args,args.region):
    for atomtype in range(1,args.types+1):
      mask = np.logical_not(np.logical_and(atoms[:,2] == atomtype, 
                                           between(atoms[:,dims[region.dim]+3],region.bounds)))
      atoms = atoms[mask]
      box[0], box[1] = region.bounds, region.bounds
  atoms[:,0] = np.arange(len(atoms))+1

  write_conf(args.output, atoms, bonds, box, types, [1] * types["atoms"], args.description)
  write_traj(args.output, np.delete(atoms[:,:6],1,1), box, mode='w')

if __name__ == '__main__':
  main()
