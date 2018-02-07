from __future__ import division
import sys
from conf_tools import *
from pbc_tools import makeBox
import numpy as np
import argparse
from numpy.random import rand
import matplotlib.pyplot as plt

MODE = 'r'

def make_size_surface(A, l, R):

  # need to calculate m and n (number of sites in the x an y directions)
  n, m = int(l[0]*np.sqrt(2/(np.sqrt(3)*A))), int(l[1]*np.sqrt(np.sqrt(3)/A/2))
  # distances between sites
  dx, dy = l[0]/n, l[1]/m
  print(dx*dy/2)
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
  surface[:,:2] -= [0.5*l[0],0.5*l[1]] #-50.0

  return surface, dx, dy

def make_density_surface(A, l, R):

  # calculate grid characteristic
  d = np.sqrt(2*A/np.sqrt(3))
  # need to calculate m and n (number of sites in the x an y directions)
  dx, dy = d, np.sqrt(3)*d/2
  n, m = int(l[0]/dx), int(l[1]/dy)
  l[0], l[1] = n*dx, m*dy
  # if there are an odd number of rows, add an extra row (m+1)
  if (m % 2):
    m -= 1
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
  surface[:,:2] -= [0.5*l[0],0.5*l[1]] #-50.0

  return surface, dx, dy

def makeFile(coords,atomtype):

  nRows = coords.shape[0]
  info = np.zeros((nRows,9))
  info[:,0] = np.arange(nRows)+1
  info[:,1] = 1
  info[:,2] = atomtype

  return info

def main():
  
  parser = argparse.ArgumentParser()
  parser.add_argument('-o','--output',default="surface")
  parser.add_argument('-r','--radius',type=float,default=0)
  parser.add_argument('-t','--atomtype',default=1,type=int)
  parser.add_argument('-f','--formats',nargs='+',choices={'lammps','conf','xyz'},
                      default='conf')
  parser.add_argument('--dims',dest='sideLengths',nargs=2,type=float,default=40)
  parser.add_argument('--packing',choices={'hexagonal','random'})
  parser.add_argument('--normal',choices={'x','y','z'},default='z')
  parser.add_argument('--density',type=float)
  parser.add_argument('--borders',action='store_true')
  parser.add_argument('--maintain',choices={'size','density'})
  args = parser.parse_args()

  a, b, c, atomRadius, nPerPart = 12.5, 12.5, 25, 1, 4684
  e = np.sqrt(1-(a/c)**2)
  if args.density is not None: areaDensity = args.density
  else: areaDensity = 2*np.pi*a**2*(1+c*np.arcsin(e)/(a*e))/nPerPart
  print(areaDensity)

  if args.packing == 'hexagonal':
    # generate a surface with a hole in the center
    if args.maintain == 'size': surface, dx, dy = make_size_surface(areaDensity, args.sideLengths, args.radius)
    elif args.maintain == 'density': surface, dx, dy = make_density_surface(areaDensity, args.sideLengths, args.radius)
    ###### VERY IMPORTANT ######
    # the rounding step is very important to getting a proper shape. 
    # open question: is this inherent to the calculation or is it an error 
    #                in the implementation?
    surface = np.around(surface,6)
    axes = {'x':0,'y':1,'z':2}
    surface[:,axes[args.normal]], surface[:,axes['z']] = surface[:,axes['z']], surface[:,axes[args.normal]].copy()
    box = makeBox(3, np.array([dx/2,dy,1.0]), surface)
    header = 'Hexagonal surface for walls'

  elif args.packing is 'random':
    box = np.array([[-22.2461,22.2461],[-22.2461,22.2461],[-1.0,1.0]])
    sides = box[:,1] - box[:,0]
    n = int(sides.prod()*0.64/np.pi*6)
    surface = sides * rand(n, 3) + box[:,0]
    info = makeFile(surface,args.atomtype)
    atoms = np.concatenate((info.astype('|S10'), surface.astype('|S10')), axis=1)
    header = 'Random surface for substrate'

  # remove sites from a the center within a circle of radius R
  surface = surface[np.sqrt(surface[:,0]**2 + surface[:,1]**2) > args.radius]
  # generate index information for configuration file
  info = makeFile(surface,args.atomtype)
  atoms = info #.astype('|S10')

  # combine coordinate and index info
  atoms[:,3:6] = surface #.astype('|S10')

  # generate input files for lammps and VMD
  if 'conf' in args.formats:
    write_conf(args.output+'_out', atoms, title=header, box=box)
  if 'xyz' in args.formats:
    write_xyz(args.output+'_out', atoms[:,2:6])
  if 'lammps' in args.formats:
    write_traj(args.output+'_out', np.delete(atoms[:,:6],1,axis=1).astype(float), box, mode='w')

  ## print out the specific surface areas for both structures
  specific_area = args.sideLengths[0]*args.sideLengths[1]/len(surface)
  x,y,z = np.diff(box,axis=1)
  specific_area = (x*y - np.pi*args.radius**2)/len(surface)
  #specific_area = np.prod(surface[:,:2].max(0)-surface[:,:2].min(0))/len(surface)
  #specific_area = np.sqrt(3)/2*(surface[:,0].max(0)-surface[:,0].min(0))**2/len(surface)
  #print('[surface]\n  a = {}\n[particle]\n  a = {}\n[error]\n  abs = {}\n  rel = {}'.format(specific_area, areaDensity, specific_area-areaDensity, (specific_area-areaDensity)/areaDensity*100))

  print('Specific Area: {}'.format(specific_area))
  print('x: {}, y: {}, z: {}'.format(x,y,z))

if __name__ == '__main__':
  main()
