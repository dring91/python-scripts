import numpy as np
from conf_tools import *
import argparse
from sys import argv,exit
from pbc_tools import PBC,unwrap
import matplotlib.pyplot as plt
from itertools import izip

def calcTransformations(frame):
  Ry = lambda a: np.array([[ np.cos(a),          0, np.sin(a)],
                           [         0,          1,         0],
                           [-np.sin(a),          0, np.cos(a)]])
  Rz = lambda a: np.array([[ np.cos(a), -np.sin(a),         0],
                           [ np.sin(a),  np.cos(a),         0],
                           [         0,          0,         1]])
  # initialize the arrays that will contain the particle centers and final rotation matrices
  rotations = np.zeros((len(frame),3,3))
  centers = np.zeros((len(frame),3))
  
  for p,particle in enumerate(frame):
    # calculate and remove the center of mass
    centers[p] = frame[p].mean(0)
    particle = frame[p] - centers[p]
    
    # determine the axes to rotate onto by finding the major and minor axes of the particle
    radii = np.sqrt(np.einsum('...i,...i',particle,particle))
    b = radii.max() # b is the major axis (ie the length)
    # xyz_b is the vector that the major axis lies along.
    xyz_b = particle[radii.argmax()]

    # Initializes a rotation as the identity matrix
    rotations[p] = np.eye(3)  
    
    # In order to rotate the particles onto a common axis, the major axis must be in the first octant
    # The first if-statement tests whether z is positive or negative and calculates a rotation of 180. 
    if xyz_b[2] < 0:
      R180 = Ry(np.pi)
      xyz_b = R180.dot(xyz_b)
      rotations[p] = rotations[p].dot(R180)
    
    # The second if-statement tests whether y is positive or negative and calculates a rotation of 180. 
    if xyz_b[1] < 0:
      R180 = Rz(np.pi)
      xyz_b = R180.dot(xyz_b)
      rotations[p] = rotations[p].dot(R180)
      
    # The third if-statement tests whether x is positive or negative and calculates a rotation of 90. 
    if xyz_b[0] < 0:
      R90 = Rz(np.pi/2)
      xyz_b = R90.T.dot(xyz_b)
      rotations[p] = R90.T.dot(rotations[p])      

    # Now, we can define the angles for the rotation matrices
    phi_b   = np.arccos(xyz_b[2]/b)
    theta_b = np.arccos(xyz_b[0]/(b*np.sin(phi_b)))
    rotations[p] = Rz(theta_b).T.dot(rotations[p])
    rotations[p] = Ry(phi_b).T.dot(rotations[p])
    
  return rotations, centers

def main():

  # compute the probability that a polymer is not in contact with two particles
  # Calculate COM of polymers
  # Compute distance between COM and each particle and generate boolean if within rc
  # Add up the number of booleans and bin into histogram of contacts
  # divide histogram by number of polymers

  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input")
  parser.add_argument("-o", "--output")
  parser.add_argument("-n", "--nFrames", type=int)
  parser.add_argument("-s", "--nSteps", type=int)
  parser.add_argument("-l", "--nMon", type=int)
  args = parser.parse_args()

  with open(args.output,"w") as out:
    out.write('')

  with open(args.input, "r") as file:
    for n in range(args.nSteps):
      for s in range(args.nFrames/args.nSteps):
        time, box, frame = readTrj(file)

      # prepare frame for analysis
      frame = frame[np.argsort(frame[:,0])]
      frame[:,2:] = frame[:,2:] * (box[:,1] - box[:,0]) + box[:,0]
      #frame[:,2:] = frame[:,2:] * (box[:,1] - box[:,0]) + box[:,0]

      # separate particle and polymer arrays and calculate box limits
      particles = frame[frame[:,1] == 3][:,2:]
      polymers = frame[frame[:,1] == 1][:,2:]
      limits = np.zeros(2)
      limits[0] = particles[:,2].min()
      limits[1] = particles[:,2].max()

      # important constants
      rc = 1.5
      a,b = 12.5 + rc,25 + rc
      nPart, nSite, nChains, nDims = 54, 4684, len(polymers)/args.nMon, 3
      
      # reshape arrays
      particles = particles.reshape((nPart, nSite, nDims))
      polymers = polymers.reshape((nChains, args.nMon, nDims))
      polymers = unwrap(polymers,box)
      com = polymers.mean(1)
      com[:,:2] = PBC(com[:,:2],com[:,:2],box[:2])

      # compute difference between chains and particles
      bridges = np.zeros(4)
      for chain in com:
        diff = particles - chain
        dist = np.sqrt(np.einsum('ijk,ijk->ij',diff,diff))
        contacts = (dist < rc).sum(1)
        touching = (contacts > 0).sum()
        bridges[touching] += 1
      probability = bridges/nChains

      with open(args.output,"a") as out:
        out.write("%d " % time)
        out.write("%f %f %f %f\n" % tuple(probability))

if __name__ == '__main__':
  main()
