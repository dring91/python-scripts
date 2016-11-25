from __future__ import division
import sys
from conf_tools import *
from pbc_tools import makeBox
import numpy as np

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

def makeFile(coords):

  nRows = coords.shape[0]
  info = np.zeros((nRows,3), dtype=int)
  info[:,0] = np.arange(nRows)+1
  info[:,1] = 1
  info[:,2] = 3

  return info

def main():
  # check for input filename from commandline
  try:
    filename = sys.argv[1]
  except IndexError:
    print "No Filename given"
    sys.exit()

  # check for cylinder radius from commandline
  try:
    R = float(sys.argv[2])
  except IndexError:
    print "Using default R value"
    R = 25.0

  a, b, c, atomRadius, nPerPart = 12.5, 12.5, 25, 1, 4684
  areaDensity = 4*np.pi*(((a*b)**1.6 + (a*c)**1.6 + (b*c)**1.6)/3)**(1/1.6)/nPerPart 
  sideLength = 4*R

  # generate a surface with a hole in the center
  surface = makeSurface(areaDensity, sideLength, R)
  ###### VERY IMPORTANT ######
  # the rounding step is very important to getting a proper shape. 
  # open question: is this inherent to the calculation or is it an error 
  #                in the implementation?
  surface = np.around(surface,6)

  box = makeBox(3, surface)

  # generate index information for configuration file
  info = makeFile(surface)
  # combine coordinate and index info
  atoms = np.concatenate((info.astype('|S10'), surface.astype('|S10')), axis=1)

  # generate input files for lammps and VMD
  write_conf(filename+'_out', atoms, title='Surface with cut-out R = %d' % R, box=box)
  write_xyz(filename+'_out', atoms[:,2:])

if __name__ == '__main__':
  main()
