from sys import argv, exit
import numpy as np
import argparse
from conf_tools import *
from pbc_tools import *

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument("-o","--output")
  parser.add_argument("-c","--nChains", type=int)
  parser.add_argument("-m","--nMon", type=int)
  parser.add_argument("--ar", type=float)
  args = parser.parse_args()

  # total number of lj beads in melt
  total = args.nMon * args.nChains

  # approximate melt density and volume
  density = 0.8
  volume = total / density
  
  # make box
  box = np.zeros((3,2))
  h = (volume / args.ar**2)**(1/3.0)
  print h
  box[2,1] = h
  box[:2,1] = h * args.ar
  box = (box.T - 0.5 * (box[:,1] + box[:,0])).T

  # 40,000 beads and 400 polymers
  # one side is approximately 40
  # xy plane is approximately 1600
  # z would be 25 planes

  bond_length = 0.97

  # perform random walks to place polymers (method from Kremer and Grest)
  # random walk with bond length = 0.97
  # restrict backfolding by: abs(r[i-1] - r[i+1]) > 1.02
  # (mind the periodic boundaries)
  # check for problems such as stretching and compute density, equilibrium R, etc.

  # Start with just one chain

  # randomly generate the bond angles of the polymer
  nAngles = args.nMon - 1
  angles = np.random.rand(args.nChains, nAngles, 2)
  angles[:,:,0] = angles[:,:,0] * (180 - 63.4) + 63.4
  angles[:,:,1] = angles[:,:,1] * 360

  # transform bond angles into polar angles
  angles[:,1:,0] = [[180 - angles[c,a,0] - angles[c,a-1,0] for a in range(1,nAngles)] for c in range(args.nChains)]
  angles[:,:,1] = np.cumsum(angles[:,:,1],axis=1)
  angles *= np.pi / 180
  
  # generate xyz coordinates
  coords = np.zeros((args.nChains,args.nMon,3))
  coords[:,1:,0] = np.sin(angles[:,:,0])*np.cos(angles[:,:,1])
  coords[:,1:,1] = np.sin(angles[:,:,0])*np.sin(angles[:,:,1])
  coords[:,1:,2] = np.cos(angles[:,:,0])
  coords *= bond_length
  coords = np.cumsum(coords, axis=0)

  # check for overlaps
  #overlap = [(np.sqrt((coords[:,c,0]-coords[:,:c,0])**2 + 
  #                    (coords[:,c,1]-coords[:,:c,1])**2 + 
  #                    (coords[:,c,2]-coords[:,:c,2])**2
  #                   ) < bond_length) for c in range(1,args.nMon)]

  # wrap the coordinates across the boundary?

  coords = coords.reshape((args.nMon * args.nChains, 3))

  # massage data into output file
  atoms = np.zeros((total,6))
  atoms[:,0] = np.arange(total)+1
  atoms[:,1] = np.repeat(np.arange(args.nChains)+1,args.nMon)
  atoms[:,2] = 1
  atoms[:,3:] = coords

  bonds = np.zeros((args.nChains * (args.nMon - 1), 4))
  bonds[:,0] = np.arange(args.nChains * (args.nMon - 1)) + 1
  bonds[:,1] = 1
  numbers = np.arange(total).reshape((args.nChains,args.nMon))
  bonds[:,2] = numbers[:,:-1].reshape((args.nChains*(args.nMon-1)))
  bonds[:,2] = numbers[:,1:].reshape((args.nChains*(args.nMon-1)))

  # output coordinates
  write_xyz(args.output, np.concatenate((np.ones((total,1)),coords), axis=1))
  write_conf(args.output, atoms, bonds, box, {"atoms":1, "bonds":1})
  write_traj(args.output, np.delete(atoms,1,axis=1), box, mode='w')
  
if __name__ == '__main__':
  main()
