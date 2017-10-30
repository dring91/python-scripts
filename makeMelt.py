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

def make_box(boundary, limits, chains, density, total, radius=None, volume=None):

  box = np.zeros((3,2))
  # Generate random starting points for chains
  start_pts = np.random.rand(chains, 3)
  volume = total / density

  if boundary == 'box':
    l = (volume)**(1/3)
    box = np.full_like(box,l/2)
    box[:,0] *= -1
    print(box)
    
    start_pts = np.random.rand(chains, 3)
    start_pts = (start_pts - 0.5) * (box[:,1] - box[:,0])
  
  if boundary == 'plane':
    # make box
    box[:2] = np.array(limits).reshape((2,2))
    h = volume / np.diff(box[:2],axis=1).prod()
    box[2,1] = h
    print(box)

    start_pts[:,:2],start_pts[:,2] = (start_pts[:,:2] - 0.5) * (box[:2,1] - box[:2,0]),\
                                     start_pts[:,2] * h
    #info = np.stack((np.arange(args.nChains)+1,np.ones(args.nChains)),axis=1)
    #write_traj('start-pts', np.concatenate((info,start_pts),axis=1), box, mode='w')

  if boundary == 'cylinder':
    # make box
    h = volume / (np.pi * radius**2)
    box[2,1] = h
    box[:2,1] = radius+0.5
    box[:2,0] = -radius-0.5

    # scale starting points to r, theta, z values
    start_pts *= [radius-0.5, 2*np.pi, h]
    # convert starting point values to x and y
    start_pts[:,0], start_pts[:,1] = start_pts[:,0]*np.cos(start_pts[:,1]), \
                                     start_pts[:,0]*np.sin(start_pts[:,1])

  return box, start_pts

def main():

  """ Algorithm
      perform random walks to place polymers (method from Kremer and Grest)
      random walk with bond length = 0.97
      restrict backfolding by: abs(r[i-1] - r[i+1]) > 1.02
  """

  # parse input
  parser = argparse.ArgumentParser()
  parser.add_argument("-o","--output",help="name of output file sans extension")
  parser.add_argument("-c","--nChains",type=int,help="Number of chains in melt")
  parser.add_argument("-m","--nMon", type=int,help="Number of monomers to a chain")
  parser.add_argument("-l","--limits",nargs='+',type=float,help="Box limits")
  parser.add_argument("-r","--radius",type=float,help="Radius of cylindrical melt")
  parser.add_argument("--reflect",default=[],nargs='+',choices={'x','y','z'})
  parser.add_argument("-b", "--boundary", choices=['plane','cylinder','box'])
  parser.add_argument("-d","--density",type=float)
  args = parser.parse_args()

  # total number of lj beads in melt
  total = args.nMon * args.nChains
  bond_length = 0.97
  dims = {'x':0,'y':1,'z':2}

  # approximate melt density and volume
  box, start_pts = make_box(args.boundary, args.limits, args.nChains, args.density, total)

  # randomly generate the bond angles of the polymer
  nAngles = args.nMon - 1
  angles = np.random.rand(args.nChains, nAngles, 2)
  bond_min = np.arccos((2*bond_length**2-1.02**2)/2*bond_length**2)
  angles[:,0,0] = angles[:,0,0] * 2 * np.pi
  angles[:,1:,0] = angles[:,1:,0] * (np.pi - bond_min) + bond_min
  angles[:,:,1] = angles[:,:,1] * 2 * np.pi

  # initialize array of cartesian coordinates
  coords = np.zeros((args.nChains, args.nMon, 3))
  # generate the bond vectors in local, relative coordinate systems
  coords[:,1:,0] = bond_length*np.sin(angles[:,:,0])*np.cos(angles[:,:,1])
  coords[:,1:,1] = bond_length*np.sin(angles[:,:,0])*np.sin(angles[:,:,1])
  coords[:,1:,2] = bond_length*np.cos(angles[:,:,0])
  # rotate vectors onto each other to obtain their local, non-relative coordinates
  coords[:,1:-1] = np.array([[rotations(rotations(head,theta,1,inv=True),
                                        np.pi-phi,2,inv=True)
                      for head, tail, (phi, theta) in zip(chain[1:-1],chain[2:],angle[1:])]
                      for chain, angle in zip(coords, angles)])

  # accumulate xyz coordinates to form a walk
  coords = np.cumsum(coords, axis=1)

  # add starting points
  coords = np.swapaxes(coords,0,1)
  coords = coords + start_pts
  coords = np.swapaxes(coords,0,1)

  # compute bond lengths
  bond_lengths = np.diff(coords,axis=1)
  bond_lengths = np.sqrt(np.sum(bond_lengths**2, axis=2))
  print(bond_lengths.min(), bond_lengths.max())

  # reshape coordinates
  coords = coords.reshape((args.nMon * args.nChains, 3))

  # deal with reflection
  types = np.ones(len(coords))
  if args.reflect != [] and args.boundary == 'cylinder':
    for _ in range(args.nMon//10):
      mask = (coords[:,0]**2+coords[:,1]**2)>(args.radius-0.5)**2
      angle = np.arctan2(coords[:,1], coords[:,0])
      boundary = (args.radius-0.5)*np.stack((np.cos(angle),np.sin(angle)),axis=1)
      coords[:,:2] -= 2*(coords[:,:2]-boundary)*np.stack((mask,mask),axis=1)
  if args.reflect != [] and args.boundary == 'plane':
    for bnd in args.reflect:
      for _ in range(args.nMon):
        above, below = coords[:,dims[bnd]] > box[dims[bnd],1], coords[:,dims[bnd]] < box[dims[bnd],0]
        coords[:,dims[bnd]] -= 2*(coords[:,dims[bnd]]-box[dims[bnd],1])*above
        coords[:,dims[bnd]] -= 2*(coords[:,dims[bnd]]-box[dims[bnd],0])*below

  # compute bond lengths
  coords = coords.reshape((args.nChains,args.nMon,3))
  bond_lengths = np.diff(coords,axis=1)
  bond_lengths = np.sqrt(np.sum(bond_lengths**2, axis=2))
  print(bond_lengths.min(), bond_lengths.max())
  coords = coords.reshape((args.nMon * args.nChains, 3))

  # massage data into output file
  atoms = np.zeros((total,9))
  atoms[:,0] = np.arange(total)+1
  #atoms[:,1] = np.repeat(np.arange(args.nChains)+1,args.nMon)
  chain = list(range(1,args.nMon//2+1)) + list(range(1,args.nMon - args.nMon//2 + 1)[::-1])
  atoms[:,1] = np.tile(chain,args.nChains)
  atoms[:,2] = types
  atoms[:,3:6] = coords

  bonds = np.zeros((args.nChains * (args.nMon - 1), 4))
  bonds[:,0] = np.arange(args.nChains * (args.nMon - 1)) + 1
  bonds[:,1] = 1
  numbers = np.arange(total).reshape((args.nChains,args.nMon)) + 1
  bonds[:,2] = numbers[:,:-1].reshape((args.nChains*(args.nMon-1)))
  bonds[:,3] = numbers[:,1:].reshape((args.nChains*(args.nMon-1)))

  diff = coords[bonds[:,2].astype(int)-1] - coords[bonds[:,3].astype(int)-1]
  lengths = np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2)

  # output coordinates
  #write_xyz(args.output, np.concatenate((np.ones((total,1)),coords), axis=1))
  write_conf(args.output, atoms, bonds, box, {"atoms":1, "bonds":1})
  #box[2] = (coords[:,2].min(), coords[:,2].max())
  write_traj(args.output, np.delete(atoms[:,:6],1,axis=1), box, mode='w')
  
if __name__ == '__main__':
  main()
