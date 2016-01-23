#!/usr/bin/python

import numpy as np
from sys import argv, exit
from getopt import *
import matplotlib.pyplot as plt

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
  rotations = {}
  centers = {}
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
  
def main():

  # tortuosity:
  #   tau = ls/l
  # 1. Construct paths through packing 
  #    - use Monte Carlo sampling (going up decreases energy)
  #    - Lay down a grid in the void space and construct paths using it
  # 2. Calculate the arc/chord ratio for all paths
  # 3. Average the ratios

  # Alternate path forming strategies:
  # - follow nearest neighbors
  # - pick endpoints and minimize distance between route and grid
  # - Bisect route and attach new node to a point in the grid within a certain distance
  # - 
  
  # Command line args processing
  trajFile, outFile, nFrames = getArgs(argv)
  
  # number of insertions
  nInsert = 1000000
  nSlices = 100
  nPaths = 200

  # define ellipsoid parameters
  r_cut = 1.625
  a = 25.0/2
  b = 50.0/2
  nAtoms = 4684
  nPart = 54
  nDim = 3

  with open(trajFile,'r') as inp:
    for n in range(nFrames):
      # read and filter data
      time, frame = readFrame(inp)
      frame = frame[frame[:,0] == 3][:,1:]

      # make box
      box = np.zeros((3,2))
      box[:,0] = frame.min(0)
      box[:,1] = frame.max(0)

      # reshape data
      frame[:,2] -= box[2,0]
      frame = frame.reshape((nPart,nAtoms,nDim))

      # unwrap particle coordinates
      frame = unwrap(frame, box)
      rotations, centers = calcTransformations(frame, box)
      
      dSlice = (box[2,1] - box[2,0]) / nSlices
      for j in range(nSlices):
        # Select insertion points
        points = np.random.rand(nInsert/nSlices,3)
        points[:,:2] = points[:,:2]*(box[:2,1] - box[:2,0]) + box[:2,0]
        points[:,2] = (points[:,2]+j)*dSlice
        
        # test for overlap
        test = np.ones_like(points[:,0],dtype=bool)
        for i,particle in enumerate(frame):
          cPoints = points - centers[i]
          rPoints = cPoints.dot(rotations[i].T)
          # particle = (particle - centers[i]).dot(rotations[i].T)
          # plt.plot(particle[:,1], particle[:,2])
          
          # test if the points satisfy the equation. If they do, they overlap a particle
          test = np.logical_and((rPoints[:,0]**2 + rPoints[:,1]**2)/(a+r_cut)**2 
                                                + rPoints[:,2]**2/(b+r_cut)**2 >= 1, test)
        # plt.scatter(points[test][:,0],points[test][:,1])
        # plt.show()
        if j == 0:
          lengths = np.zeros(nPaths) 
          chords = np.zeros((nPaths,3))
          paths = np.zeros((nPaths,nSlices,3))
          indices = np.random.choice(len(points[test]),nPaths)
          paths[:,0] = points[test][indices]
        else:
          for p,path in enumerate(paths):
            # Calculate distance between path and all points
            displacement = points[test]-path[j-1]
            distance = np.linalg.norm(displacement,axis=1)
            # pick minimum distance and add to path length
            lengths[p] += np.abs(distance.min())
            # add displacement to path
            chords[p] += displacement[np.argmin(distance)]
            paths[p,j] = points[test][np.argmin(distance)]
            
      tau = lengths / np.linalg.norm(chords,axis=1)
      print 2*tau.mean()**2
      # plt.plot(tau)
      # plt.show()
       
      ######## Checking ######## 
      # plt.plot(dSlice*(np.arange(nSlices) + 0.5) / (box[2,1] - box[2,0]),paths)
      ax1 = plt.subplot(231)
      ax2 = plt.subplot(234)
      ax3 = plt.subplot(132)
      ax4 = plt.subplot(133)
      ax1.plot(2*tau**2)
      for p in range(nPaths):
        # plot 3 views (xy, xz, yz)
        # plt.plot(np.linalg.norm(paths[p,:,:2],axis=1),paths[p,:,2])
        ax2.plot(paths[p,:,0],paths[p,:,1])
        ax2.set(aspect='equal')
        ax2.set_xlim((-50,50))
        ax2.set_ylim((-50,50))
        ax3.plot(paths[p,:,0],paths[p,:,2])
        ax3.set_xlim((-50,50))
        ax4.plot(paths[p,:,1],paths[p,:,2])
        ax4.set_xlim((-50,50))
        # print paths[p,:,2]
      # plt.xlabel('Height, z (LJ)')
      # plt.ylabel(r'porosity, $\phi$ (v/v)')
      # plt.title(r'$\phi_{bulk}$ = %.2f, $\phi_{avg}$ = %.2f' % (porosity_bulk, porosity_avg))
      plt.tight_layout()
      plt.show()
      ##########################    

    # ####### Output ########
    # # zSlices = [dSlice*(k + 0.5) / (box[2,1] - box[2,0]) for k in range(nSlices)]
    # with open(outFile, 'w') as otp:
    #   for path in paths:
    #     otp.write('#   x    y    z\n')
    #     [otp.write('%f %f %f\n' % tuple(l)) for l in path]
    #     otp.write('\n')
    # #######################
      
    print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
