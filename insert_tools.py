import numpy as np
from sys import argv, exit
from future_builtins import zip

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

def calc_porosity(frame, box, ns, params):
  # initialize porosity
  porosity = np.zeros((ns['nBins'],2))
  porosity[:,0] = (np.arange(ns['nBins']) + 0.5) / ns['nBins']

  # calculate bin size
  frame = frame[:,1:].reshape((ns['nPart'],ns['nAtoms'],ns['nDim']))
  dSlice = (box[2,1] - box[2,0]) / ns['nBins']

  # unwrap particle coordinates
  frame = unwrap(frame, box)

  # find the transforms to translate and rotate particles
  rotations, centers = calcTransformations(frame, box)
  
  # Select insertion points
  points = np.random.rand(ns['nInsert'],3)
  points = points*(box[:,1] - box[:,0]) + box[:,0]
  
  # test for overlap
  test = np.zeros_like(points[:,0],dtype=bool)
  for c,r in zip(centers, rotations):
    test = insert(points, c, r, box, test, params)

  # partition points by vertical height
  bins = ((points[:,2] - box[2,0]) / dSlice).astype(int)
  porosity[:,1] = np.zeros(ns['nBins'])

  # calculate the ratio of accepted to total insertions
  for i in range(ns[3]):
    insertions = np.logical_not(test[bins == i])
    porosity[i,1] = insertions.sum() / float(len(insertions))

  return porosity

def main():
  
  # Command line args processing
  args = getArgs(argv)
  
  outFile = '%s_N5_%d_%d' % (args.output, args.particles, args.bins)

  params = {'a':25*0.5, 'b':50*0.5, 'cut-off':0}
  ns = {'nPart':54, 'nAtoms':4684, 'nDim':3, 'nBins':args.bins, 'nInsert':args.insertions}

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

      porosity = calc_porosity(frame, box, ns, params)
      
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
