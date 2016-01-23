#!/bin/python

import numpy as np
from sys import argv, exit
from getopt import *
import matplotlib.pyplot as plt
# from numpy.stats import mode

def getArgs(argv):
  
  # defaults
  trajFile = 'Ellipsoids/613766.ocoee.nics.utk.edu/N_5_undersat_traj.xyz'
  outFile = ''
  steps = 0

  try:
    opts, args = getopt(argv[1:],'t:o:n:')
  except GetoptError:
    print 'randomWalk.py -t <filename> -o <filename> -n <steps>'
    exit(2)
  for opt, arg in opts:
    if opt == '-t':
      trajFile = arg
    elif opt == '-o':
      outFile = arg
    elif opt == '-n':
      steps = arg

  return trajFile, outFile, int(steps)
  
def readFrame(file):
  line = file.readline().strip()
  nAtoms = int(line)
  atoms = np.zeros((nAtoms,4))

  line = file.readline()
  time = int(line.split()[-1])

  for i in range(nAtoms):
    line = file.readline().split()
    atoms[i,:] = np.array(line, dtype=float) 

  return time, atoms
   
def PBC(dx,x,interval):
  mask = (dx < interval[:,0])
  x = x*mask + 2*interval[:,1]*mask + x*np.logical_not(mask)
  mask = (dx >= interval[:,1])
  x = x*mask + 2*interval[:,0]*mask + x*np.logical_not(mask)
    
  return x
  
def unwrap(coords, box):
  for i,p in enumerate(coords):
    coords[i,:,:2] = PBC(p[:,:2]-p[0,:2],p[:,:2],box[:2,:])
      
  return coords
  
def calcTransformations(frame, box):
  # define rotation matrices
  Ry = lambda a: np.array([[ np.cos(a),          0, np.sin(a)],
                           [         0,          1,         0],
                           [-np.sin(a),          0, np.cos(a)]])
  Rz = lambda a: np.array([[ np.cos(a), -np.sin(a),         0],
                           [ np.sin(a),  np.cos(a),         0],
                           [         0,          0,         1]])
  rotations = np.zeros((len(frame),3,3))
  centers = np.zeros((len(frame),3))
  for p,particle in enumerate(frame):
    # calculate and remove the center of mass
    COM = frame[p].mean(0)
    centers[p] = COM
    particle = frame[p] - COM
    
    # determine the axes to rotate onto
    radii = np.sqrt(np.einsum('...i,...i',particle,particle))
    b = radii.max()
    xyz_b = particle[radii.argmax()]

    rotations[p] = np.eye(3)  
    
    if xyz_b[2] < 0:
      R180 = Ry(np.pi)
      xyz_b = R180.dot(xyz_b)
      rotations[p] = rotations[p].dot(R180)
    
    if xyz_b[1] < 0:
      R180 = Rz(np.pi)
      xyz_b = R180.dot(xyz_b)
      rotations[p] = rotations[p].dot(R180)
      
    if xyz_b[0] < 0:
      R90 = Rz(np.pi/2)
      xyz_b = R90.T.dot(xyz_b)
      rotations[p] = R90.T.dot(rotations[p])      

    # define angles
    phi_b   = np.arccos(xyz_b[2]/b)
    theta_b = np.arccos(xyz_b[0]/(b*np.sin(phi_b)))
    rotations[p] = Rz(theta_b).T.dot(rotations[p])
    rotations[p] = Ry(phi_b).T.dot(rotations[p])
    
  return rotations, centers
  
def writeHeader(filename, header, mode):
  with open(filename, mode) as file:
    file.write(header)

def writeFile(filename, data, time, mode):
  with open(filename, mode) as file:
    file.write('# time: %d\n' % time)
    file.write('#   z (sigma)   porosity\n')
    [file.write('%f %f\n' % l) for l in data]
    file.write('\n')

def writeLine(filename, line):
  with open(filename, 'a') as file:
    file.write('%f %f\n' % line)

def transformCheck(frame, centers, rotations):
  for i,p in enumerate(frame):
    # plt.plot( (p-centers[i]).dot(rotations[i])[:,0], 
    #           (p-centers[i]).dot(rotations[i])[:,2] )
    f, axes = plt.subplots()
    plt.plot( rotations[i].dot((p-centers[i]).T).T[:,0], 
              rotations[i].dot((p-centers[i]).T).T[:,2] )
    axes.set(aspect='equal')
  plt.show()

def main():
  
  # Command line args processing
  trajFile, outFile, nFrames = getArgs(argv)
  
  # number of insertions
  nInsert = 1000000
  outFile = '%s_N5_%d' % (outFile, nInsert)

  # define ellipsoid parameters
  max_r = 40
  sizes = 0.5*np.arange(1,max_r)
  a = 25.0/2
  b = 50.0/2
  nAtoms = 4684
  nPart = 54
  nDim = 3
  
  header = '# Porosity for N = 5 oversaturated infiltration\n'
  # header += '# 12 frames calculated\n'
  writeHeader(outFile, header, 'w')
  mode = 'a'

  with open(trajFile,'r') as inp:
    # read and filter data
    time, frame = readFrame(inp)
    frame = frame[frame[:,0] >= 3]

  # make box
  box = np.zeros((3,2))
  box[:,0] = frame[:,1:].min(0)
  box[:,1] = frame[:,1:].max(0)
  Lz = box[2,1] - box[2,0]
  ends = np.array([0.20, -0.20]) * Lz
  box[2,:] += ends
  print box

  # reshape data
  frame = frame[:,1:].reshape((nPart,nAtoms,nDim))

  # unwrap particle coordinates
  frame = unwrap(frame, box)

  # find the transforms to translate and rotate particles
  rotations, centers = calcTransformations(frame, box)
  print 'Transforms calculated'
  # transformCheck(frame, centers, rotations)

  for r_cut in sizes:
    # randomly select points for insertion
    points = np.random.rand(nInsert,3)
    points = points*(box[2,1] - box[2,0]) + box[:,0]
    
    # test for overlap
    accepted = 0
    test = np.zeros_like(points[:,0],dtype=bool)
    for i,particle in enumerate(frame):
      cPoints = points - centers[i]
      cPoints[:,:2] = PBC(cPoints[:,:2], cPoints[:,:2], box[:2,:])
      rPoints = cPoints.dot(rotations[i].T)
      
      # test if the points satisfy the equation. If they do, they overlap a particle
      test = np.logical_or((rPoints[:,0]**2 + rPoints[:,1]**2)/(a+r_cut)**2 
                                            + rPoints[:,2]**2 /(b+r_cut)**2 < 1, test)

    # calculate the ratio of accepted to total insertions
    accepted = sum(np.logical_not(test))
    print '%d insertions accepted' % accepted
    accepted = accepted / float(nInsert)
    
    writeLine(outFile, (r_cut, accepted))
    print 'Size %f written to file' % r_cut
    
  print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
