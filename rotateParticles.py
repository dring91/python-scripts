def calcTransformations(frame):
  # input: particle xyz data
  # output: rotation matrices to rotate the nanoparticles onto a common axis, centers of mass
  # purpose: This function specifically calculates rotation matrices and angles that specify the locations of 
  #          ellipsoidal nanoparticles relative to a common axis. The function is used for insertion so that
  #          test particles inside of the ellipsoids can be cleanly rejected. It may be possible to use the
  #          rotation angles to define the quaternion for a particle.

  # define rotation matrices
  # (the subscript is the axis of rotation)
  Ry = lambda a: np.array([[ np.cos(a),          0, np.sin(a)],
                           [         0,          1,         0],
                           [-np.sin(a),          0, np.cos(a)]])
  Rz = lambda a: np.array([[ np.cos(a), -np.sin(a),         0],
                           [ np.sin(a),  np.cos(a),         0],
                           [         0,          0,         1]])
  # initialize the arrays that will contain the particle centers and final rotation matrices
  rotations = np.zeros((len(frame),3,3))
  centers = np.zeros((len(frame),3))
  
  # The shape of every frame is (nParticles, nBeads, nDimensions)
  # Therefore, the first loop is over the number of nanoparticles, 
  #     the second over the number of beads per nanoparticle, etc.
  #
  # schematically,
  #
  #      particle 1    ...   particle M
  #    
  # [ [[ x1, y1, z1 ],     [ x1, y1, z1 ],
  #    [ x2, y2, z2 ],     [ x2, y2, z2 ],
  #    ...                 ...
  #    [ xN, yN, zN ]] ... [ xN, yN, zN ]] ] 

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
