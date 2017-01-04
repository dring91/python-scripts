import numpy as np
from conf_tools import *
import argparse
from sys import argv,exit
from pbc_tools import PBC,unwrap
import matplotlib.pyplot as plt

def ellipsoid(coords, axes=[1,1,1]):

  A = np.eye(3)*np.array(axes,dtype=float)**(-2)

  A = np.einsum('ij,ki->kj',A,coords)
  A = np.einsum('ij,ij',A,coords)

  return A

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

  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input")
  parser.add_argument("-o", "--output")
  parser.add_argument("-n", "--nFrames", type=int)
  parser.add_argument("-s", "--nSteps", type=int)
  args = parser.parse_args()

  with open(args.output,"w") as out:
    out.write('')

  with open(args.input, "r") as file:
    for n in range(args.nSteps):
      for s in range(args.nFrames/args.nSteps):
        time, box, frame = readTrj(file)

      frame = frame[np.argsort(frame[:,0])]
      frame[:,2:] = frame[:,2:] * (box[:,1] - box[:,0]) + box[:,0]
      #frame[:,2:] = frame[:,2:] * (box[:,1] - box[:,0]) + box[:,0]
      particles = frame[frame[:,1] == 3][:,2:]
      polymers = frame[frame[:,1] == 1][:,2:]
      limits = np.zeros(2)
      limits[0] = particles[:,2].min()
      limits[1] = particles[:,2].max()

      # reshape particle array
      particles = particles.reshape((54,4684,3))
      upper30 = (particles[:,:,2] > (limits[1]-30)).sum()
      lower30 = (particles[:,:,2] < (limits[0]+30)).sum()
      middle = np.logical_and(particles[:,:,2] >= (limits[0]+30),particles[:,:,2] <= (limits[1]-30)).sum()

      rc = 1.5
      a,b = 12.5 + rc,25 + rc
      nPart, nSite = 54, 4684
      
      # select only monomers within a given shell
      not_close = [0,0,0,0]
      particles = unwrap(particles,box)
      rotations, centers = calcTransformations(particles)
      for p in range(nPart):
        diff = polymers - centers[p]
        diff[:,:2] = PBC(diff[:,:2],diff[:,:2],box[:2])
        diff = diff.dot(rotations[p].T)
        shell = polymers[(diff[:,0]**2+diff[:,1]**2)/a**2+diff[:,2]**2/b**2 <= 1]
        shell[:,:2] = PBC(shell[:,:2]-centers[p,:2],shell[:,:2],box[:2])
        #plt.plot(particles[p,:,0],particles[p,:,1],'o')
        #plt.plot(shell[:,0],shell[:,1],'o')
        #plt.show()
        for site in particles[p]:
          diff = site - shell
          dist = np.sqrt(np.einsum('...i,...i',diff,diff))
          zeros = (0 == (dist <= rc).sum())
          not_close[0] += zeros
          if site[2] > (limits[1]-30): not_close[1] += zeros
          if site[2] < (limits[0]+30): not_close[2] += zeros
          if np.logical_and(site[2] <= (limits[1]-30),site[2] >= (limits[0]+30)): 
            not_close[3] += zeros
      probability = [not_close[0] / float(nPart * nSite), 
                     not_close[1] / float(upper30), 
                     not_close[2] / float(lower30),
                     not_close[3] / float(middle)]

      with open(args.output,"a") as out:
        out.write("%d %f %f %f %f\n" % tuple([time]+probability))

if __name__ == '__main__':
  main()
