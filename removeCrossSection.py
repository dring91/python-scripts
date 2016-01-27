#!/bin/python

import numpy as np
from sys import argv, exit
from getopt import *

NDIMS = 3
AREA = 100*100

def getArgs(argv):
  
  try:
    opts, args = getopt(argv[1:],'t:o:n:i:s:')
  except GetoptError:
    print """randomWalk.py -t <filename> -o <filename> -n <steps> 
                           -i <nInsert> -s <nSlices>"""
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

  return trajFile, outFile, int(steps), int(insertions), int(slices)
  
def readFrame(file,nCols=4):
  line = file.readline().strip()
  nAtoms = int(line)

  line = file.readline()
  time = int(line.split()[-1])

  atoms = np.fromfile(file, float, nAtoms*nCols, ' ')
  atoms = atoms.reshape((nAtoms,nCols))

  return time, atoms

def interpolate(x, low, high):
  if x == low[0]:
    return low[1]
  elif x == high[0]:
    return high[1]
  elif x < high[0] and x > low[0]:
    slope = (high[1] - low[1])/(high[0] - low[0])
    return slope*(x - low[0]) + low[1]
  else:
    print 'x out of range'
    exit(0)

def integrateHeight(data, stopIntegrationAt):
  height = data[0,0]
  prev = data[0]
  for row in data[1:]:
    if row[1] > stopIntegrationAt >= prev[1]:
      height = interpolate(stopIntegrationAt, prev[::-1], row[::-1])
      break
    prev = row
  
  return height
 
def makeHistogram(atoms,nBins,limits):
  binSize = (limits[1]-limits[0])/float(nBins)
  hist = np.zeros((nBins,2))
  hist[:,0] = (np.arange(nBins)+0.5)*binSize / limits[1]
  atoms = atoms[np.logical_and(atoms >= limits[0], atoms < limits[1])]
  bins = ((atoms - limits[0]) / nBins).astype(int)
  # loop through bins and add atoms
  for bin in range(nBins):
    hist[bin,1] = (bins == bin).sum()
    # try:
    # except IndexError:
    #   print "Particle outside box"
  hist[:,1] /= float(AREA*binSize)
  
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
    particle = frame[p] - centers[p]
    
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

def calcPorosity(particles, box, nInsert, nSlices, params):
  porosity = np.zeros((nSlices,2))
  porosity[:,0] = (np.arange(nSlices) + 0.5) / nSlices
  dSlice = (box[2,1] - box[2,0]) / nSlices
  particles = particles.reshape((params[4],params[3],NDIMS))

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
    test = insert(points, centers[i], rotations[i], box, test, params[:3])

  # partition points by vertical height
  bins = ((points[:,2] - box[2,0]) / dSlice).astype(int)

  # calculate the ratio of accepted to total insertions
  for j in range(nSlices):
    insertions = np.logical_not(test[bins == j])
    porosity[j,1] = insertions.sum() / float(len(insertions))

  return porosity

def findHeight(particles, box, nInsert, nSlices, params, polymers, stopIntegrationAt):
  porosity = calcPorosity(particles, box, nInsert, nSlices, params)
  # How does this work with findHeight?
  density = makeHistogram(polymers, nSlices, box[2])
  density[:,1] /= porosity[:,1]
  # Integrate density within packing
  height = integrateHeight(density, stopIntegrationAt)

  return density, height

def main():
  
  # Command line args processing
  trajFile, outFile, nFrames, nInsert, nSlices = getArgs(argv)
  
  outFile = '%s_%d' % (outFile, nInsert)

  # params = (a,b,r_cut,nAtomsPerPart,nPart)
  params = (25*0.5, 50*0.5, 0, 4684, 54)

  with open('density_'+outFile, 'w') as otp:
    otp.write('# Density profiles\n')
  with open('height_'+outFile, 'w') as otp:
    otp.write('# height profiles at two cut-offs\n')

  with open(trajFile,'r') as inp:
    for n in range(nFrames):
      # read data and separate by atomtypes
      time, frame = readFrame(inp)
      particles = frame[frame[:,0] >= 3][:,1:]
      polymers = frame[frame[:,0] == 1][:,3]
      pmin = polymers.min()

      # make box
      box = np.zeros((NDIMS,2))
      box[:,0] = particles.min(0)
      box[:,1] = particles.max(0)

      # Calculate important system values
      # remove all scaling; rescale at the end if necessary
      lowerBound = pmin - box[2,0]
      nOutside = (polymers < box[2,0]).sum()
      totalDensity = (polymers < box[2,1]).sum() / (box[2,1] - pmin) / AREA
      bulkDensity = nOutside / (-lowerBound) / AREA
      
      cutOff85 = 0.85
      cutOff99 = 0.99
      density85 = None
      density99 = None
      if (-bulkDensity) * lowerBound >= cutOff99 * totalDensity * (box[2,1] - lowerBound):
        height85 = totalDensity / bulkDensity * cutOff85 * (box[2,1] - lowerBound) + lowerBound
        height99 = totalDensity / bulkDensity * cutOff99 + lowerBound
      elif bulkDensity * lowerBound >= cutOff85 * totalDensity * (box[2,1]+lowerBound):
        height85 = totalDensity / bulkDensity * cutOff85 + lowerBound
        stopIntegrationAt = totalDensity * cutOff99 + bulkDensity * lowerBound
        density99, height99 = findHeight(particles, box, nInsert, nSlices, params, polymers, stopIntegrationAt)
      else:
        stopIntegrationAt = totalDensity * cutOff85 + bulkDensity * lowerBound
        density85, height85 = findHeight(particles, box, nInsert, nSlices, params, polymers, stopIntegrationAt)
        stopIntegrationAt = totalDensity * cutOff99 + bulkDensity * lowerBound
        density99, height99 = findHeight(particles, box, nInsert, nSlices, params, polymers, stopIntegrationAt)

      density = np.ones((nSlices,3))
      if density99 is None:
        density = [[0, bulkDensity, bulkDensity]]
      elif density85 is None:
        density[:,0] = density99[:,0]
        density[:,1] *= bulkDensity
        density[:,2] = density99[:,1]
      else:
        density[:,:2] = density85
        density[:,2] = density99[:,1]

      # output density and height values
      with open('density_'+outFile, 'a') as otp:
        header = '# time: %d\n#  z  density85  density99' % time
        np.savetxt(otp, density, fmt='%.5f', header=header, footer='\n', comments='')
      with open('height_'+outFile, 'a') as otp:
        np.savetxt(otp, [[time, height85, height99]], fmt='%.5f',comments='')
      
  print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()

# Count number of atoms outside the packing
# - if the fraction of atoms outside the packing exceeds the cut-off,
#   then use method 1 to calculate height
# - if the fraction is less, then use method 2
# 
# Method 1: height = rho_T/rho_b*C + lower
# Method 2: Adjust cut-off with rho_T*C - rho_b*(-lower) and integrate rho
#           within the packing
