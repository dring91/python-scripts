from sys import argv, exit
import numpy as np
import argparse
from conf_tools import *
from pbc_tools import *
import matplotlib.pyplot as plt

def rotations(v, angle, axis, inv=False):

  axes = np.arange(3)
  [i,j] = axes[axes != axis]
  if inv: i,j = j,i

  # what is the best way to perform broadcasting?
  R = np.zeros((3,3))
  R[axis,axis] = 1
  R[i,i] = np.cos(angle) 
  R[j,j] = np.cos(angle) 
  R[i,j] = np.sin(angle) 
  R[j,i] = -np.sin(angle) 

  return v.dot(R.T)

def main():

  ## Algorithm
  # perform random walks to place polymers (method from Kremer and Grest)
  # random walk with bond length = 0.97
  # restrict backfolding by: abs(r[i-1] - r[i+1]) > 1.02

  # parse input
  parser = argparse.ArgumentParser()
  parser.add_argument("-o","--output")
  parser.add_argument("-c","--nChains", type=int)
  parser.add_argument("-m","--nMon", type=int)
  parser.add_argument("--ar", type=float)
  parser.add_argument("-r","--radius",type=float)
  parser.add_argument("--reflect",action='store_true')
  parser.add_argument("-b", "--boundary", choices=['box','cylinder'])
  args = parser.parse_args()

  # total number of lj beads in melt
  total = args.nMon * args.nChains
  bond_length = 0.97

  # approximate melt density and volume
  density = 0.8
  volume = total / density
  
  if args.boundary == 'box':
    # make box
    box = np.zeros((3,2))
    h = (volume / args.ar**2)**(1/3.0)
    box[2,1] = h
    box[:2,1] = h * args.ar
    box = (box.T - 0.5 * (box[:,1] + box[:,0])).T

    # Generate random starting points for chains
    start_pts = np.random.rand(args.nChains, 3)
    start_pts = (start_pts - 0.5) * (box[:,1] - box[:,0])

  if args.boundary == 'cylinder':
    # make box
    box = np.zeros((3,2))
    h = volume / (np.pi * args.radius**2)
    box[2,1] = h
    box[:2,1] = args.radius
    box[:2,0] = -args.radius

    # Generate random starting points for chains
    start_pts = np.random.rand(args.nChains, 3)
    start_pts *= [args.radius-0.5, 2*np.pi, h]
    start_pts[:,0], start_pts[:,1] = start_pts[:,0]*np.cos(start_pts[:,1]), \
                                     start_pts[:,0]*np.sin(start_pts[:,1])

  # randomly generate the bond angles of the polymer
  nAngles = args.nMon - 1
  angles = np.random.rand(args.nChains, nAngles, 2)
  bond_min = np.arccos((2*bond_length**2-1.02**2)/2*bond_length**2)
  angles[:,:,0] = angles[:,:,0] * (np.pi - bond_min) + bond_min
  angles[:,:,1] = angles[:,:,1] * 2 * np.pi

  # initialize array of cartesian coordinates
  coords = np.zeros((args.nChains, args.nMon, 3))
  # generate the first bond vector
  coords[:,1,0] = bond_length*np.sin(angles[:,0,0])*np.cos(angles[:,0,1])
  coords[:,1,1] = bond_length*np.sin(angles[:,0,0])*np.sin(angles[:,0,1])
  coords[:,1,2] = bond_length*np.cos(angles[:,0,0])
  # rotate vectors onto each other to obtain their cartesian coordinates
  for i in range(args.nChains):
    for j in range(2,args.nMon):
      coords[i,j] = rotations(
                    rotations(coords[i,j-1],angles[i,j-1,1],1,inv=True),
                                            np.pi-angles[i,j-1,0],2,inv=True)
 
  # accumulate xyz coordinates to form a walk
  coords = np.cumsum(coords, axis=1)

  # add starting points
  coords = np.swapaxes(coords,0,1)
  coords = coords + start_pts
  coords = np.swapaxes(coords,0,1)

  # deal with reflection
  if args.reflect:
    for _ in range(args.nMon):
      mask = (coords[:,:,0]**2+coords[:,:,1]**2)>(args.radius-0.5)**2
      coords -= 2*(coords-args.radius+0.5)*np.stack((mask,mask,mask),axis=2)

  # reshape coordinates
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
  numbers = np.arange(total).reshape((args.nChains,args.nMon)) + 1
  bonds[:,2] = numbers[:,:-1].reshape((args.nChains*(args.nMon-1)))
  bonds[:,3] = numbers[:,1:].reshape((args.nChains*(args.nMon-1)))

  # output coordinates
  #write_xyz(args.output, np.concatenate((np.ones((total,1)),coords), axis=1))
  #write_conf(args.output, atoms, bonds, box, {"atoms":1, "bonds":1})
  write_traj(args.output, np.delete(atoms,1,axis=1), box, mode='w')
  
if __name__ == '__main__':
  main()
