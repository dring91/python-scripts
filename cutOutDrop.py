#!/usr/bin/python

from sys import argv, exit
from getopt import *
import numpy as np

def write_xyz(filename, atoms, nChains, chainLength):
  with open(filename+'.xyz','w') as otp:
    otp.writelines([str(len(atoms)),'\n','Atoms\n'])
    for line in atoms:
      otp.write('1 ')
      xyz = [str(x)+' ' for x in line]
      otp.writelines(xyz)
      otp.write('\n')

def write_conf(filename,
               atoms,
               nChains,
               title='Entangled polymer simulation with chains of N = 250\n',
               chainLength=250,
               types=[1,1], # [atomtypes,bondtypes]
               box=[[1,1],[1,1],[1,1]],
               masses=[1]):
  
  ##########################################  
  # Header comment line:
  # 
  # ###### atoms
  # ###### bonds
  # 
  # # atom types
  # # bond types
  #
  # # # xlo xhi
  # # # ylo yhi
  # # # zlo zhi
  # 
  # Masses
  #
  # # #
  # 
  # Atoms
  # 
  # atom_id molecule_id atomtype x y z
  # 
  # Bonds
  #
  # bond_id bondtype atom1 atom2
  #
  ###########################################
  if types[1] > 0:
    bonds = (chainLength-1)*nChains
  else:
    bonds = 0
  
  with open(filename+'.conf','w') as otp:  
    # write title line
    otp.write(title)
    
    # skip line
    otp.write('\n')
    
    # write number of atoms and number of bonds
    otp.write(str(len(atoms))+' atoms\n')
    if types[1] > 0:
      otp.write(str(bonds)+' bonds\n')

    # skip line
    otp.write('\n')
    
    # write types
    otp.write(str(types[0])+' atom types\n')
    if types[1] > 0:
      otp.write(str(types[1])+' bond types\n')
    
    # skip line
    otp.write('\n')
    
    # write box dimensions
    dim = 'xyz'
    [otp.write('%f %f %slo %shi\n' % (box[l][0],box[l][1],dim[l],dim[l])) for l in range(len(box))]    

     # skip line
    otp.write('\n')

    # write masses
    otp.write('Masses\n')
    otp.write('\n')
    [otp.write('%d %d\n' % (i+1,masses[i])) for i in range(len(masses))]

    # skip line    
    otp.write('\n')
    
    # write atoms
    otp.write('Atoms\n')
    otp.write('\n')
    shift = np.median(range(1,chainLength+1))
    num = lambda x: shift-np.abs((x+chainLength-1)%chainLength+1-shift)
    sym = False
    if sym:
      [otp.write('%d %d %d %f %f %f\n' % (i+1,
                                          num(i+1),
                                          types[0]-1,
                                          atoms[i][0],atoms[i][1],atoms[i][2])
                                          ) 
                                          for i in range(nChains*chainLength)]
      [otp.write('%d %d %d %f %f %f\n' % (i+1+nChains*chainLength,
                                          chainLength/2+1,
                                          types[0],
                                          atoms[i][0],atoms[i][1],atoms[i][2])
                                          )
                                          for i in range(len(atoms)-nChains*chainLength)]
    else:
      [otp.write('%d %d %d %f %f %f\n' % (i+1,
                                          np.floor(i/chainLength)+1,
                                          types[0],
                                          atoms[i][0],atoms[i][1],atoms[i][2])
                                          ) 
                                          for i in range(nChains*chainLength)]
      """
      [otp.write('%d %d %d %f %f %f\n' % (i+1+nChains*chainLength,
                                          nChains+1,
                                          types[0],
                                          atoms[i][0],atoms[i][1],atoms[i][2])
                                          )
                                          for i in range(len(atoms)-nChains*chainLength)]
                                          """

    if types[1] > 0:
      # skip line    
      otp.write('\n')
      
      # write bonds
      otp.write('Bonds\n')
      otp.write('\n')
      [[otp.write('%d %d %d %d\n' % (j*(chainLength-1)+i+1,
                                     types[1],
                                     j*(chainLength)+i+1,
                                     j*(chainLength)+i+2)
                                     )
                                     for i in range(chainLength-1)]
                                     for j in range(nChains)]

def inRange(point, interval=[0,1], left=True, right=False):
  
  if point > interval[0] and point < interval[1]:
    return True
  elif left and point == interval[0]:
    return True
  elif right and point == interval[1]:
    return True
  else:
    return False
      
def inBox(point, box=[[0,1],[0,1],[0,1]], left=True, right=False):
  
  for i in range(len(point)):
    if inRange(point[i],box[i]) == False:
      return False
  
  return True

def inSphere(point, R=[0,1], inner=True, outer=True, hemisphere=False):
  
  r = np.sqrt(np.dot(point,point))
  if hemisphere and inRange(r,R,inner,outer) and point[2] > 0:
    return True
  elif inRange(r,R,inner,outer):
    return True
  
  return False

def inCylinder(point, R, L, origin=[0,0,0], outer=True, inner=True):

  R = [0.0,R]
  L = [-0.5*L,0.5*L]
  point = np.array(point)
  origin = np.array(origin)
  r = np.sqrt(np.dot(point[1:]-origin[1:],point[1:]-origin[1:]))
  if inRange(r,R,inner,outer) and inRange(point[0]-origin[0],L):
    return True

  return False

def PBC(dx,x,interval):
  if dx <= interval[0]:
    x = x + 2*interval[1]
  elif dx > interval[1]:
    x = x + 2*interval[0]

  return x  

def calcCOM(atoms):

  COM = 0
  for atom in atoms:
    COM += atom
  COM /= len(atoms)

  return COM

def getArgs(argv):

  try:
    opts, args = getopt(argv[1:],'i:o:l:r:z:')
  except GetoptError:
    print 'randomWalk.py -i <file> -o <file> -l <chainLength> -r <radius> -z <z-origin>'
    exit(2)
  for opt, arg in opts:
    if opt == '-i':
      inFile = arg
    elif opt == '-o':
      outFile = arg
    elif opt == '-l':
      chainLength = arg
    elif opt == '-r':
      radius = arg
    if opt == '-z':
      origin = arg
 
  return inFile, outFile, int(chainLength), float(radius), float(origin)

def boundAtoms(coords, dim):
  bounds = [coords[0][dim],coords[0][dim]]
  for coord in coords:
    if coord[dim] < bounds[0]:
      bounds[0] = coord[dim]
    elif coord[dim] > bounds[1]:
      bounds[1] = coord[dim]

  return bounds

def readConf(file, atype):

  atoms = []
  box = []
  bonds = []
  header = 'header'
  for line in file:
    L = line.split()
    if len(L) > 0 and L[0] in set(['Atoms','Bonds']):
      header = L[0]
    if len(L) > 0 and L[-1] in set(['xhi','yhi','zhi']):
      box.append(L[:2])
    elif len(L) > 2 and L[2] in set(atype) and header == 'Atoms':
      atoms.append(L)
    elif len(L) > 2 and header == 'Bonds':
      bonds.append(L)
      
  return box, atoms, bonds

def cut(atoms, R, origin, chainLength):
  n = 0
  molecule = []
  drop = []
  # loop over each molecule
  for atom in atoms:
    if n == chainLength:
      # unnecessary because of unwrap
      # molecule = PBC_box(molecule, box)
      # calculate it's center of mass
      molecule = np.array(molecule)
      COM = [calcCOM(molecule[:,d]) for d in range(len(molecule[0]))]
      # check whether it's in the region
      if inRange(COM[0],[-0.5*R,0.5*R]): # inCylinder(COM, R, 0.5*R, origin):
        # add it to drop if it is
        drop.extend(molecule)
      n = 0
      molecule = []
    n += 1
    molecule.append(atom)

  return drop

def main():

  inFile, outFile, chainLength, R, z = getArgs(argv)

  try:
    with open(inFile+'.conf', 'r') as inp:
      box, atoms, bonds = readConf(inp,['1'])
  except IOError:
    print 'Cannot read input file'
    exit(0)

  atoms = [[float(i) for i in atom[3:]] for atom in atoms]
  drop = cut(atoms, R, [0,0,z], chainLength)

  # make new box
  box = [boundAtoms(drop,d) for d in range(3)]

  write_xyz(outFile, drop, len(drop)/chainLength, chainLength)
  # write_conf(filename,atoms,bonds,title,types,box,masses):
  write_conf(outFile, 
             drop,
             len(drop)/chainLength,
             'undersaturated infiltration simulation\n', 
             chainLength,
             [1,1], 
             box,
             [1]
            )

if __name__ == "__main__":
  main()
