from sys import argv, exit
import numpy as np
import argparse
from conf_tools import *
from pbc_tools import *
import matplotlib.pyplot as plt

def rotations(v, angle, axis, inv=False):

  """ 
  Routine to rotate a vector by constructing the appropriate 3D rotation matrix
      v: the vector to rotate
  angle: the angle to rotate v by
   axis: the axis to rotate v around (0:x,1:y,2:z)
    inv: switch to determine whether the rotation array is an inverse rotation
  """

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
  parser.add_argument("-r","--radius", type=float, help='only valid with cylindrical geom')
  parser.add_argument("--ar", type=float)
  parser.add_argument("--files",nargs='+',choices={'xyz','conf','lammps'})
  parser.add_argument("--geometry",choices={'cartesian','cylindrical','spherical'})
  parser.add_argument("--reflect",action="store_true",
                      help="determine behavior at boundaries")
  args = parser.parse_args()

  # total number of lj beads in melt
  total = args.nMon * args.nChains

  # approximate melt density and volume
  density = 0.8
  volume = total / density
  bond_length = 0.97

  if args.geometry == 'cartesian':
  
    # make box
    h = (volume / args.ar**2)**(1/3.0)
    box = np.zeros((3,2))
    box[2,1] = h
    box[:2,1] = h * args.ar
    box = (box.T - 0.5 * (box[:,1] + box[:,0])).T

    # Generate random starting points for chains
    start_pts = np.random.rand(args.nChains, 3)
    start_pts = (start_pts - 0.5) * (box[:,1] - box[:,0])

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
      # the third element is the first coordinate specified in a reltive axis
      for j in range(2,args.nMon):
        coords[i,j] = rotations(
                      rotations(coords[i,j-1],angles[i,j-1,1],1,inv=True),
                                np.pi-angles[i,j-1,0],2,inv=True)
 
  if args.geometry == 'cylindrical':

    # make box
    h = volume / (np.pi*(args.radius-0.5)**2)
    box = np.zeros((3,2))
    box[2,1] = h
    box[:2,1] = args.radius-0.5
    box[:2,0] = 0.5-args.radius

    # Generate random starting points for chains
    start_pts = np.random.rand(args.nChains, 3)
    start_pts *= [args.radius-0.5, 2*np.pi, h]
    start_pts[:,0], start_pts[:,1] = start_pts[:,0]*np.cos(start_pts[:,1]), \
                                     start_pts[:,0]*np.cos(start_pts[:,1])
    # calculate the minimum bond angle
    bond_min = np.arccos((2*bond_length**2-1.02**2)/(2*bond_length**2))
    # Generate angles randomly
    angles = np.random.rand(args.nMon - 1, args.nChains, 2)
    angles[:,:,0] = angles[:,:,0] * (np.pi - bond_min) + bond_min
    angles[:,:,1] = angles[:,:,1] * 2 * np.pi
    # Convert angles to coordinates
    coords = np.zeros((args.nMon, args.nChains, 3))
    coords[1,:,0] = bond_length*np.sin(angles[0,:,0])*np.cos(angles[0,:,1])
    coords[1,:,1] = bond_length*np.sin(angles[0,:,0])*np.sin(angles[0,:,1])
    coords[1,:,2] = bond_length*np.cos(angles[0,:,0])

    for i in range(2,args.nMon):
      for j in range(args.nChains):
        coords[i,j] = rotations(
                      rotations(coords[i-1,j],angles[i-1,j,1],1,inv=True),
                                np.pi-angles[i-1,j,0],2,inv=True)
      # Reflect monomers at boundaries
      if args.reflect: 
        coords[i,:,1:] -= 2*coords[i,:,1:]*(coords[i,:,0]**2 + coords[i,:,1]**2 > (args.radius-0.5)**2)
 
  # add start points
  coords = np.cumsum(coords, axis=0)
  coords += start_pts
  coords = np.swapaxes(coords,0,1)
  bonds = np.diff(coords,axis=0)
  print(np.sqrt(bonds[:,:,0]**2+bonds[:,:,1]**2+bonds[:,:,2]**2))
  coords = coords.reshape((args.nMon * args.nChains, 3))

  # massage data into output file
  atoms = np.zeros((total,9))
  atoms[:,0] = np.arange(total)+1
  atoms[:,1] = np.repeat(np.arange(args.nChains)+1,args.nMon)
  atoms[:,2] = 1
  atoms[:,3:6] = coords

  bonds = np.zeros((args.nChains * (args.nMon - 1), 4))
  bonds[:,0] = np.arange(args.nChains * (args.nMon - 1)) + 1
  bonds[:,1] = 1
  numbers = np.arange(total).reshape((args.nChains,args.nMon)) + 1
  bonds[:,2] = numbers[:,:-1].reshape((args.nChains*(args.nMon-1)))
  bonds[:,3] = numbers[:,1:].reshape((args.nChains*(args.nMon-1)))

  # output coordinates
  if 'xyz' in args.files: 
    write_xyz(args.output, np.concatenate((np.ones((total,1)),coords), axis=1))
  if 'conf' in args.files: 
    write_conf(args.output, atoms, bonds, box, {"atoms":1, "bonds":1})
  if 'lammps' in args.files: 
    write_traj(args.output, np.delete(atoms[:,:6],1,axis=1), box, mode='w')
  
if __name__ == '__main__':
  main()
