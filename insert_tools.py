import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit

from pbc_tools import *

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
    centers[p] = particle.mean(0)
    particle = particle - centers[p]
    
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
    
    ##### Transform check ####
    #f, axes = plt.subplots()
    #plt.plot( particle[:,0], 
    #          particle[:,2], '.' )
    #plt.plot( particle.dot(rotations[p].T)[:,0], 
    #          particle.dot(rotations[p].T)[:,2], '.' )
    #axes.set(aspect='equal')
    #plt.show()
    ##########################

  return rotations, centers
  
def insert(insertions, center, rotation, box, accepted, par):
  # DO NOT CHANGE TO insertions -= center
  insertions = insertions - center
  # apply minimum image convention to points before rotation
  insertions[:,:2] = PBC(insertions[:,:2], insertions[:,:2], box[:2,:])
  insertions = insertions.dot(rotation.T)

  # #### point rotation check ####
  # f, axes = plt.subplots()
  # plt.plot( insertions[:,0], insertions[:,2], '.')
  # # plt.plot( (particle-center).dot(rotation.T)[:,0],
  # #           (particle-center).dot(rotation.T)[:,2], '.')
  # axes.set(aspect='equal')
  # plt.show()
  # ##############################
  
  # test if the points satisfy the equation. If they do, they overlap a particle
  par['r'] = 0
  accepted = np.logical_or((insertions[:,0]**2 
                          + insertions[:,1]**2)/(par['a']+par['r'])**2 
                          + insertions[:,2]**2/(par['b']+par['r'])**2 < 1, accepted)

  return accepted

if __name__ == '__main__':
  main()
