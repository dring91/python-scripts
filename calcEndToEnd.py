import numpy as np
from sys import argv, exit
from argparse import ArgumentParser
from conf_tools import iTrj
import matplotlib.pyplot as plt

def main():

  # Command line args processing
  parser = ArgumentParser()
  parser.add_argument('-i', '--input')
  parser.add_argument('-o', '--output')
  parser.add_argument('-l', '--length', type=int)
  parser.add_argument('-t', '--types', type=int, nargs='+')

  args = parser.parse_args()
  
  with open(args.input,'r') as inp:
    for (time, nAtoms, box, atoms) in iTrj(inp):
      ## sort and filter atoms
      atoms = atoms[np.argsort(atoms[:,0])]
      polymers = atoms[atoms[:,1] == args.types[0]]
      cylinder = atoms[atoms[:,1] == args.types[1]]

      ## split atoms into regions to compute
      bounds = [box[2,0], cylinder[4].min(), box[2,1]]
      regions = {'inCyl':polymers[polymers[:,4] < bounds[1]], 
                 'inFilm':polymers[polymers[:,4] >= bounds[1]]}

      ## compute the following for each region
      for region in regions: # modify for dictionary iteration
        ## reshape array
        nMon = len(polymers)
        polymers = polymers.reshape((nMon / args.length, args.length, 5))

        ## calculate end-to-end distance
        differences = polymers[:,1:,2:5] - polymers[:,:-1,2:5]
        end_to_end = differences.sum(1)
        R = np.sqrt(np.einsum('ij,ij->i',end_to_end,end_to_end))
        # final array has dimensions of num_polymers

        ## Calculate bond lengths
        lengths = np.sqrt(np.einsum('ijk,ijk->ij',differences,differences))
        lengths = lengths.reshape(args.length - 1) * nMon / args.length)
        # final array has dimensions of num_bonds

        ## calculate radius of gyration
        polymers = polymers.swapaxes(0,1)
        center_of_mass = polymers[:,2:5].mean(0)
        radius_gyration = polymers[:,2:5] - center_of_mass
        radius_gyration = np.einsum('ijk,ijk->ij',radius_gyration,radius_gyration)
        radius_gyration = np.sqrt(radius_gyration.mean(0))
        # final array has dimensions of num_polymers
    
        ## concatenate the arrays to write out
        output = np.concatenate((center_of_mass, end_to_end, radius_gyration))

        with open(args.output, 'a') as otp:
        #  np.savetxt(otp)
     
  print('Finished analyzing trajectory')

if __name__ == '__main__':
  main()
