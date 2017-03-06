from __future__ import division
import sys
from conf_tools import *
from pbc_tools import makeBox
import numpy as np
import argparse
from numpy.random import rand
import matplotlib.pyplot as plt

MODE = 'r'

def makeSurface(A, l, R):

  # need to calculate m and n (number of sites in the x an y directions)
  n, m = int(l*np.sqrt(2/(np.sqrt(3)*A))), int(l*np.sqrt(np.sqrt(3)/A/2))
  # distances between sites
  dx, dy = l/n, l/m
  # if there are an odd number of rows, add an extra row (m+1)
  if (m % 2):
    m += 1
  row = np.arange(n) * dx
  row_offset = row + 0.5*dx
  surface = np.zeros((m*n,3))
  for i in range(m):
    if i % 2 == 0:
      surface[i*n:(i+1)*n,0] = row
      surface[i*n:(i+1)*n,1] = i*dy
    else:
      surface[i*n:(i+1)*n,0] = row_offset
      surface[i*n:(i+1)*n,1] = i*dy

  # shift center to position (0,0)
  surface[:,:2] -= 0.5*l #-50.0
  # remove sites from a the center within a circle of radius R
  surface = surface[np.sqrt(surface[:,0]**2 + surface[:,1]**2) > R]

  return surface

def makeFile(coords,atomtype):

  nRows = coords.shape[0]
  info = np.zeros((nRows,9), dtype=int)
  info[:,0] = np.arange(nRows)+1
  info[:,1] = 1
  info[:,2] = atomtype

  return info

def main():
  
  parser = argparse.ArgumentParser()
  parser.add_argument('-o','--output',default="surface")
  parser.add_argument('-r','--radius',type=float)
  parser.add_argument('-t','--atomtype',default=3,type=int)
  parser.add_argument('-f','--formats',nargs='+',choices={'lammps','conf','xyz'},
                      default='conf')
  args = parser.parse_args()

  a, b, c, atomRadius, nPerPart = 12.5, 12.5, 25, 1, 4684
  areaDensity = 4*np.pi*(((a*b)**1.6 + (a*c)**1.6 + (b*c)**1.6)/3)**(1/1.6)/nPerPart 
  sideLength = 40.0 # 4*args.radius

  if args.radius is not None:
    # generate a surface with a hole in the center
    surface = makeSurface(areaDensity, sideLength, args.radius)
    ###### VERY IMPORTANT ######
    # the rounding step is very important to getting a proper shape. 
    # open question: is this inherent to the calculation or is it an error 
    #                in the implementation?
    surface = np.around(surface,6)

    box = makeBox(3, np.array([0.4,0.93,0]), surface)
    # box = np.zeros((3,2))
    # box[:2] = np.array([-20,20])

  else:
    box = np.array([[-22.2461,22.2461],[-22.2461,22.2461],[-1.0,1.0]])
    sides = box[:,1] - box[:,0]
    n = int(sides.prod()*0.64/np.pi*6)
    surface = sides * rand(n, 3) + box[:,0]

    info = makeFile(surface,args.atomtype)
    atoms = np.concatenate((info.astype('|S10'), surface.astype('|S10')), axis=1)

  # generate index information for configuration file
  info = makeFile(surface,args.atomtype)
  atoms = info.astype('|S10')
  # combine coordinate and index info
  atoms[:,3:6] = surface.astype('|S10')

  # generate input files for lammps and VMD
  if 'conf' in args.formats:
    write_conf(args.output+'_out', atoms, title='Random surface for substrate', box=box)
  if 'xyz' in args.formats:
    write_xyz(args.output+'_out', atoms[:,2:6])
  if 'lammps' in args.formats:
    write_traj(args.output+'_out', np.delete(atoms[:,:6],1,axis=1).astype(float), box, mode='w')

if __name__ == '__main__':
  main()
