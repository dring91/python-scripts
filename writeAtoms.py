#!/usr/bin/python

def write_xyz(filename, atoms):
  with open(filename+'.xyz','w') as otp:
    otp.write('%d\nAtoms\n' % len(atoms))
    try:
      [otp.write('%s %s %s %s\n' % tuple(line)) for line in atoms]
    except TypeError:
      [otp.write('%d %f %f %f\n' % tuple(line)) for line in atoms]

def write_conf(filename,
               atoms,
               bonds,
               title='Entangled polymer simulation with chains of N = 250\n',
               types=[1,1], # [atomtypes,bondtypes]
               box=[[1,1],[1,1],[1,1]],
               masses=[1]):
  
  ##########################################  
  # Header comment line:
  # 
  # ###### atoms
  # ###### bonds
  # 
  # # atom types
  # # bond types
  #
  # # # xlo xhi
  # # # ylo yhi
  # # # zlo zhi
  # 
  # Masses
  #
  # # #
  # 
  # Atoms
  # 
  # atom_id molecule_id atomtype x y z
  # 
  # Bonds
  #
  # bond_id bondtype atom1 atom2
  #
  ###########################################

  with open(filename+'.conf','w') as otp:  
    # write title line
    otp.write(title+'\n')
    
    # skip line
    otp.write('\n')
    
    # write number of atoms and number of bonds
    otp.write(str(len(atoms))+' atoms\n')
    if types[1] > 0:
      otp.write(str(len(bonds))+' bonds\n')

    # skip line
    otp.write('\n')
    
    # write types
    otp.write(str(types[0])+' atom types\n')
    if types[1] > 0:
      otp.write(str(types[1])+' bond types\n')
    
    # skip line
    otp.write('\n')
    
    # write box dimensions
    dim = 'xyz'
    try:
      [otp.write('%s %s %slo %shi\n' % (box[l][0],box[l][1],dim[l],dim[l])) for l in range(len(box))]    
    except TypeError:
      [otp.write('%f %f %slo %shi\n' % (box[l][0],box[l][1],dim[l],dim[l])) for l in range(len(box))]    

     # skip line
    otp.write('\n')

    # write masses
    otp.write('Masses\n')
    otp.write('\n')
    try:
      [otp.write('%s %s\n' % (i+1,masses[i])) for i in range(len(masses))]
    except TypeError:
      [otp.write('%d %d\n' % (i+1,masses[i])) for i in range(len(masses))]

    # skip line    
    otp.write('\n')
    
    # write atoms
    otp.write('Atoms\n')
    otp.write('\n')
    # shift = np.median(range(1,chainLength+1))
    # num = lambda x: shift-np.abs((x+chainLength-1)%chainLength+1-shift)
    try:
      [otp.write(' '.join(line)+'\n') for line in atoms]
    except TypeError:
      [otp.write('%d %d %d %f %f %f\n' % tuple(line)) for line in atoms]

    if types[1] > 0:
      # skip line    
      otp.write('\n')
      
      # write bonds
      otp.write('Bonds\n')
      otp.write('\n')
      try:
        [otp.write(' '.join(line)+'\n') for line in bonds]
      except TypeError:
        [otp.write('%d %d %d %d\n' % tuple(line)) for line in bonds]