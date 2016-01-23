#!/bin/python

import numpy as np
from sys import argv, exit
from getopt import *
import matplotlib.pyplot as plt
# from numpy.stats import mode

def getArgs(argv):
  
  try:
    opts, args = getopt(argv[1:],'t:o:n:i:s:d:')
  except GetoptError:
    print 'randomWalk.py -t <filename> -o <filename> -n <steps> 
                         -i <nInsert> -s <nSlices> -d <binWidth>'
    exit(2)
  for opt, arg in opts:
    if opt == '-t':
      trajFile = arg
    elif opt == '-o':
      outFile = arg
    elif opt == '-n':
      steps = arg
    elif opt == '-i':
      insertions = arg
    elif opt == '-s':
      slices = arg
    elif opt == '-d':
      delta = arg

  return trajFile, outFile, int(steps), int(insertions), int(slices), float(delta)
  
def readFrame(file,nCols=4):
  line = file.readline().strip()
  nAtoms = int(line)

  line = file.readline()
  time = int(line.split()[-1])

  atoms = np.fromfile(file, float, nAtoms*nCols, ' ')
  atoms = atoms.reshape((nAtoms,nCols))

  return time, atoms

def makeHistogram(atoms,nBins,lower):
  area = 100*100
  hist = np.zeros((nBins,2))
  binSize = (1-lower)/float(nBins)
  hist[:,0] = (np.arange(nBins)+0.5)*binSize+lower
  atoms = atoms[np.logical_and(atoms > lower, atoms <= 1)
  
  for atom in atoms:
    if lower < atom <= 1:
      bin = int((atom-lower)/binSize)
      try:
        hist[bin,1] += 1
      except IndexError:
        print "Particle outside box"
  hist[:,1] /= area*binSize
  
  return hist
  
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
    centers[p] = frame[p].mean(0)
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
    
    # #### Transform check ####
    #   # plt.plot( (p-centers[i]).dot(rotations[i])[:,0], 
    #   #           (p-centers[i]).dot(rotations[i])[:,2] )
    #   f, axes = plt.subplots()
    #   plt.plot( rotations[i].dot((p-centers[i]).T).T[:,0], 
    #             rotations[i].dot((p-centers[i]).T).T[:,2] )
    #   axes.set(aspect='equal')
    # plt.show()
    # #########################

  return rotations, centers

def insert(insertions, center, rotation, box, accepted, (a,b,r)):
  # DO NOT CHANGE TO insertions -= center
  insertions = insertions - center
  # apply minimum image convention to points before rotation
  insertions[:,:2] = PBC(insertions[:,:2], insertions[:,:2], box[:2,:])
  insertions = insertions.dot(rotation.T)

  # test if the points satisfy the equation. If they do, they overlap a particle
  accepted = np.logical_or((insertions[:,0]**2 + insertions[:,1]**2)/(a+r)**2 
                                               + insertions[:,2]**2/(b+r)**2 < 1, accepted)

  return accepted

def calcPorosity(particles, box, nInsert, nSlices, params, nAtoms, nPart, nDim):
  dSlice = (box[2,1] - box[2,0]) / nSlices
  particles = particles.reshape((nPart,nAtoms,nDim))

  # unwrap particle coordinates
  particles = unwrap(particles, box)

  # find the transforms to translate and rotate particles
  rotations, centers = calcTransformations(particles, box)
  
  # Select insertion points
  points = np.random.rand(nInsert,3)
  points = points*(box[:,1] - box[:,0]) + box[:,0]
  
  # test for overlap
  test = np.zeros_like(points[:,0],dtype=bool)
  for i,particle in enumerate(particles):
    test = insert(points, centers, rotations, box, test, params)

  # partition points by vertical height
  bins = ((points[:,2] - box[2,0]) / dSlice).astype(int)
  porosity[:,1] = np.zeros(nSlices)

  # calculate the ratio of accepted to total insertions
  for j in range(nSlices):
    accepted = np.logical_not(test[bins == j])
    porosity[j,1] = insertions.sum() / float(len(insertions))

  return porosity

def main():
  
  # Command line args processing
  trajFile, outFile, nFrames, nInsert, nSlices, delta = getArgs(argv)
  
  porosity = np.zeros((nSlices,2))
  porosity[:,0] = (np.arange(nSlices) + 0.5) / nSlices
  outFile = '%s_N5_%d' % (outFile, nInsert)

  # params = (a,b,r_cut)
  params = (25*0.5, 50*0.5, 0)
  nAtoms = 4684
  nPart = 54
  nDim = 3

  with open(outFile, 'w') as otp:
    otp.write('# Porosity for N = 5 oversaturated infiltration\n')

  with open(trajFile,'r') as inp:
    for n in range(nFrames):
      # read and filter data
      time, frame = readFrame(inp)
      particles = frame[frame[:,0] >= 3][:,1:]

      # make box
      box = np.zeros((3,2))
      box[:,0] = particles.min(0)
      box[:,1] = particles.max(0)

      porosity = calcPorosity(particles, box, nInsert, nSlices, params, nAtoms, nPart, nDim)
 
      polymers = frame[frame[:,0] == 1][:,3]

      polymers -= box[2,0]
      polymers /= (box[2,1] - box[2,0])
      lower = (polymers.min() - box[2,0])/(box[2,1] - box[2,0])
      density = makeHistogram(polymers, nSlices, lower)
      density[:,1] /= (box[2,1] - box[2,0])
      
      with open(outFile, 'a') as otp:
        np.savetxt(otp, 
                   ,
                   fmt='%.5f',
                   header='# time: %d\n#   z (sigma)   porosity\n' % time,
                   footer='\n',
                   comments='')
      
    print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()

"""
n = n_b + n_p
n = rho_b * v_b + rho_p * v_p
rho = rho_b * (h_b/L) + rho_p * (h_p/L)

Histgram diagram

Top of packing
-----------------------------------------|8 
                                         |7
                                         |6
                                         |5
                                         |4
                                         |3
                                         |2
Bottom of packing                        |1
-----------------------------------------|0 
Top of Polymer bulk                      |1 
                                         |2 t dependent
                                         |3
-----------------------------------------|4
Bottom of Polymer bulk

packing doesn't change height, but the polymer bulk does
Is it necessary to bin below the bottom of the packing?
"""
