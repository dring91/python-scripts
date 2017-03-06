from __future__ import division
import sys
from conf_tools import *
from pbc_tools import makeBox
import numpy as np
import argparse

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
  info = np.zeros((nRows,9), dtype=int)
  info[:,0] = np.arange(nRows)+1
  info[:,1] = 1
  info[:,2] = 3

  return info

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('-o','--output',default="cylinder")
  parser.add_argument('-r','--radius',type=float)
  parser.add_argument('-l','--length',type=float,default=100)
  parser.add_argument('-t','--atomtype',default=3,type=int)
  parser.add_argument('-f','--formats',nargs='+',choices={'lammps','conf','xyz'},
                      default='conf')
  args = parser.parse_args()

  a, b, c, atomRadius, nPerPart = 12.5, 12.5, 25, 1, 4684
  areaDensity = 4*np.pi*(((a*b)**1.6 + (a*c)**1.6 + (b*c)**1.6)/3)**(1/1.6)/nPerPart 

  # generate a cylinder with radius R and length L
  cylinder = makeCylinder(areaDensity, (args.radius,args.length))
  ###### VERY IMPORTANT ######
  # the rounding step is very important to getting a proper shape. 
  # open question: is this inherent to the calculation or is it an error 
  #                in the implementation?
  cylinder = np.around(cylinder,6)

  box = makeBox(3, np.array([1.0,1.0,0]), cylinder)

  # generate indices, molecule, and atom types
  info = makeFile(cylinder)
  atoms = info.astype('|S10')
  # combine coordinate and index information
  atoms[:,3:6] = cylinder.astype('|S10')

  # generate input files for VMD and LAMMPS
  if 'conf' in args.formats:
    write_conf(args.output+'_out', atoms, title='Capillary R = %d' % args.radius, box=box)
  if 'lammps' in args.formats:
    write_traj(args.output+'_out', np.delete(atoms[:,:6],1,axis=1).astype(float), box, mode='w')
  if 'xyz' in args.formats:
    write_xyz(args.output+'_out', atoms[:,2:])

if __name__ == '__main__':
  main()
