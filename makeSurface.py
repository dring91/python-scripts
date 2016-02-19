#!/usr/bin/python

from __future__ import division
import sys
from conf_tools import *
from pbc_tools import makeBox
import numpy as np
from math import pi

MODE = 'r'

def makeSurface(A, l):

  # need to calculate m and n
  n, m = int(l*np.sqrt(2/(np.sqrt(3)*A))), int(l*np.sqrt(np.sqrt(3)/A/2))
  dx, dy = l/n, l/m
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

  return surface

def makeFile(coords):

  nRows = coords.shape[0]
  info = np.zeros((nRows,3), dtype=int)
  info[:,0] = np.arange(nRows)+1
  info[:,1] = 1
  info[:,2] = 3

  return info

def main():
  try:
    filename = sys.argv[1]
  except IndexError:
    print "No Filename given"
    sys.exit()

  try:
    R = float(sys.argv[2])
  except IndexError:
    print "Using default R value"
    R = 25

  atype = ['1','2']
  nMon, nPerPart = 50, 4684
  a, b, c, atomRadius = 12.5, 12.5, 25, 1
  areaDensity = 4*np.pi*(((a*b)**1.6 + (a*c)**1.6 + (b*c)**1.6)/3)**(1/1.6)/nPerPart 
  sideLength = 100

  surface = makeSurface(areaDensity, sideLength)
  surface = np.around(surface,6)
  surface[:,:2] += -50.0
  surface = surface[np.sqrt(surface[:,0]**2 + surface[:,1]**2) > R]

  box = makeBox(3, surface)

  info = makeFile(surface)
  atoms = np.concatenate((info.astype('|S10'), surface.astype('|S10')), axis=1)

  write_conf(filename+'_out', atoms, title='Surface with cut-out R = %d' % R, box=box)
  write_xyz(filename+'_out', atoms[:,2:])

if __name__ == '__main__':
  main()
