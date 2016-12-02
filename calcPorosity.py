import numpy as np
from sys import argv, exit
from future_builtins import zip
import argparse

from conf_tools import readTrj
from pbc_tools import *

def getArgs(argv):

  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input")
  parser.add_argument("-o", "--output")
  parser.add_argument("-f", "--frames", type=int)
  parser.add_argument("-p", "--particles", type=int)
  parser.add_argument("-b", "--bins", type=int)
  
  return parser.parse_args()
  
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
    
    # #### Transform check ####
    # f, axes = plt.subplots()
    # plt.plot( particle.dot(rotations[p].T)[:,0], 
    #           particle.dot(rotations[p].T)[:,2] )
    # axes.set(aspect='equal')
    # plt.show()
    # #########################

  return rotations, centers
  
def insert(insertions, center, rotation, box, accepted, (a,b,r)):
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
  accepted = np.logical_or((insertions[:,0]**2 + insertions[:,1]**2)/(a+r)**2 
                                               + insertions[:,2]**2/(b+r)**2 < 1, accepted)

  return accepted

def main():
  
  # Command line args processing
  #trajFile, outFile, nFrames, nInsert, nSlices = getArgs(argv)
  args = getArgs(argv)
  
  porosity = np.zeros((args.bins,2))
  porosity[:,0] = (np.arange(args.bins) + 0.5) / args.bins
  outFile = '%s_N5_%d_%d' % (args.output, args.particles, args.bins)

  # params = (a,b,r_cut)
  params = (25*0.5, 50*0.5, 0)
  nAtoms = 4684
  nPart = 54
  nDim = 3

  with open(args.output, 'w') as otp:
    otp.write('# Porosity for N = 5 oversaturated infiltration\n')

  with open(args.input,'r') as inp:
    for n in range(args.frames):
      # read and filter data from xyz file
      #time, frame = readFrame(inp)
      #frame = frame[frame[:,0] >= 3]

      # read and filter data from traj file
      time, box, frame = readTrj(inp)
      frame = frame[frame[:,1] == 3]

      # make box
      box = np.zeros((3,2))
      box[:,0] = frame[:,2:].min(0)
      box[:,1] = frame[:,2:].max(0)
      dSlice = (box[2,1] - box[2,0]) / args.bins

      frame = frame[:,1:].reshape((nPart,nAtoms,nDim))

      # unwrap particle coordinates
      frame = unwrap(frame, box)

      # find the transforms to translate and rotate particles
      rotations, centers = calcTransformations(frame, box)
      
      # Select insertion points
      points = np.random.rand(args.particles,3)
      points = points*(box[:,1] - box[:,0]) + box[:,0]
      
      # test for overlap
      test = np.zeros_like(points[:,0],dtype=bool)
      for c,r in zip(centers, rotations):
        test = insert(points, c, r, box, test, params)

      # partition points by vertical height
      bins = ((points[:,2] - box[2,0]) / dSlice).astype(int)
      porosity[:,1] = np.zeros(nSlices)

      # calculate the ratio of accepted to total insertions
      for i in range(nSlices):
        insertions = np.logical_not(test[bins == i])
        porosity[i,1] = insertions.sum() / float(len(insertions))
       
      with open(outFile, 'a') as otp:
        np.savetxt(otp, 
                   porosity,
                   fmt='%.5f',
                   header='# time: %d\n#   z (sigma)   porosity\n' % time,
                   footer='\n',
                   comments='')
      
    print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
