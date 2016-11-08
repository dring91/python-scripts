#!/usr/bin/python

from __future__ import division
import sys
from conf_tools import *
from pbc_tools import makeBox
import numpy as np

MODE = 'r'

def ring(r,n,period):
  dth = 2*np.pi/n
  coords = np.zeros((n,3))
  for i in range(n):
    coords[i,:2] = [r*np.cos(dth*(i+period)), r*np.sin(dth*(i+period))]

  return coords

def makeCylinder(A, (r,h)):

  # need to calculate m and n (the number of sites in the axial and polar directions)
  m, n = int(h*np.sqrt(2/(np.sqrt(3)*A))), int(np.pi*r*np.sqrt(2*np.sqrt(3)/A))
  # L is the distance between sites in the polar direction
  L = 2*np.pi*r/n
  # generate two rings that will form each repeated section of the cylinder
  layer = ring(r,n,0)
  layer_offset = ring(r,n,0.5)
  # initialize the cylinder
  cylinder = np.zeros((m*n,3))
  # loop over the number of rings
  for i in range(m):
    if i % 2 == 0:
      cylinder[i*n:(i+1)*n,:2] = layer[:,:2]
      cylinder[i*n:(i+1)*n,2] = i*L*np.sqrt(3)/2
    else:
      cylinder[i*n:(i+1)*n,:2] = layer_offset[:,:2]
      cylinder[i*n:(i+1)*n,2] = i*L*np.sqrt(3)/2

  return cylinder

def makeFile(coords):

  nRows = coords.shape[0]
  info = np.zeros((nRows,3), dtype=int)
  info[:,0] = np.arange(nRows)+1
  info[:,1] = 1
  info[:,2] = 3

  return info

def main():
  # get output filename from commandline
  try:
    filename = sys.argv[1]
  except IndexError:
    print "No Filename given"
    sys.exit()

  a, b, c, atomRadius, nPerPart = 12.5, 12.5, 25, 1, 4684
  areaDensity = 4*np.pi*(((a*b)**1.6 + (a*c)**1.6 + (b*c)**1.6)/3)**(1/1.6)/nPerPart 

  # get cylinder dimensions from commandline
  try:
    cylDims = (float(sys.argv[2]),100)
  except IndexError:
    cylDims = (25,100)
    print "Default R value used"

  # generate a cylinder with radius R and length L
  cylinder = makeCylinder(areaDensity, cylDims)
  ###### VERY IMPORTANT ######
  # the rounding step is very important to getting a proper shape. 
  # open question: is this inherent to the calculation or is it an error 
  #                in the implementation?
  cylinder = np.around(cylinder,6)

  box = makeBox(3, cylinder)

  # generate indices, molecule, and atom types
  info = makeFile(cylinder)
  # combine coordinate and index information
  atoms = np.concatenate((info.astype('|S10'), cylinder.astype('|S10')), axis=1)

  # generate input files for VMD and LAMMPS
  write_conf(filename+'_out', atoms, title='Capillary R = %d' % cylDims[0], box=box)
  write_xyz(filename+'_out', atoms[:,2:])

if __name__ == '__main__':
  main()
