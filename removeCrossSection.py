import numpy as np
from sys import argv, exit
from getopt import *
from conf_tools import readTrj

NDIMS = 3
AREA = 100*100
RC = 1

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

def integrateHeight(data, cutOff):
  # rescale and accumulate data values
  data[1,1] *= -2*data[0,0]
  data[2:,1] *= (data[2,0] - data[1,0])
  data[:,1] = np.cumsum(data[:,1]) / data[:,1].sum()
  above = data[data[:,1] > cutOff]
  below = data[data[:,1] <= cutOff]
 
  if below[-1,1] == cutOff:
    height = below[-1,0]
  else:
    height = interpolate(cutOff, below[-1,::-1], above[0,::-1])
  
  return height
 
def makeHistogram(atoms,nBins,limits):
  # make a histogram with different sized bins
  hist = np.zeros((nBins+2,2))
  bottomSize = limits[0] - atoms.min()
  binSize = (limits[1]-limits[0])/float(nBins)
  hist[0,0] = -bottomSize # limits[0] - 0.5 * bottomSize
  hist[1:,0] = np.arange(nBins+1)*binSize # + limits[0]
  hist[1,1] = (atoms < limits[0]).sum()

  atoms = atoms[np.logical_and(atoms >= limits[0], atoms < limits[1])]
  bins = ((atoms - limits[0]) / binSize).astype(int)
  # loop through bins and add atoms
  for bin in range(nBins-2):
    hist[bin+2,1] = (bins == bin+2).sum()
  hist[:,1] /= float(AREA)
  hist[2:,1] /= binSize
  hist[1,1] /= bottomSize
  
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
  porosity[:,0] = np.arange(nSlices) / nSlices
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

def findHeight(particles, box, nInsert, nSlices, params, polymers, cutOff):
  # Calculate the packing porosity
  porosity = calcPorosity(particles, box, nInsert, nSlices, params)
  # calculate the superficial density of the polymers in the packing
  density = makeHistogram(polymers, nSlices, box[2])
  # divide by the porosity to recovery the real density
  density[2:,1] /= porosity[:,1]
  # Integrate real density to the appropriate cut-off
  height = integrateHeight(np.copy(density), cutOff)

  return density, height, porosity

def main():
  
  # Command line args processing
  trajFile, outFile, nFrames, nInsert, nSlices = getArgs(argv)
  
  outFile = '%s_%d' % (outFile, nInsert)

  # params = (a,b,r_cut,nAtomsPerPart,nPart)
  params = (25*0.5, 50*0.5, RC, 4684, 54)

  with open('density_'+outFile, 'w') as otp:
    otp.write('# Density profiles\n')
  with open('height_'+outFile, 'w') as otp:
    otp.write('# height profiles at two cut-offs\n')

  with open(trajFile,'r') as inp:
    for n in range(nFrames):
      # read data and separate by atomtypes
      time, box, frame = readTrj(inp)
      frame = frame[frame[:,0].argsort()]
      particles = frame[frame[:,1] >= 3][:,2:]
      polymers = frame[frame[:,1] == 1][:,4]
      pmin = polymers.min()

      cutOff = 0.85
      density, height85, porosity = findHeight(particles, box, nInsert, nSlices, params, polymers, cutOff)
      cutOff = 0.99
      height99 = integrateHeight(np.copy(density), cutOff)

      # rescale height
      dmax = density[-1,0]
      density[:,0] /= dmax
      height85 /= dmax
      height99 /= dmax
      # concatenate porosity and density???

      # output density and height values
      with open('density_'+outFile, 'a') as otp:
        header = '# time: %d\n#  z  density  porosity' % time
        np.savetxt(otp, density, fmt='%.5f', header=header, footer='\n', comments='')
      with open('height_'+outFile, 'a') as otp:
        np.savetxt(otp, [[time, height85, height99, density[0,0]]], fmt='%d %.5f %.5f %.5f',comments='')
      
  print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
