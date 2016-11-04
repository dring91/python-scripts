import numpy as np

def readConf(file, atype):

  atoms = []
  box = []
  bonds = []
  header = 'header'
  for line in file:
    L = line.split()
    if len(L) > 0 and L[0] in ['Atoms','Bonds']:
      header = L[0]
    if len(L) > 0 and L[-1] in ['xhi','yhi','zhi']:
      box.append(L[:2])
    elif len(L) > 2 and L[2] in atype and header == 'Atoms':
      atoms.append(L)
    elif len(L) > 2 and header == 'Bonds':
      bonds.append(L)
      
  return box, atoms, bonds

def readTrj(file, nCols=5):

  """ default """
  time = None
  nAtoms = 0
  box = None 
  atoms = None
  item = None

  line = file.readline().strip()
  while (line != ''):
    L = line.split()
    if L[0] == 'ITEM:':
      item = L[1]
    if item is not None:
      if item == 'TIMESTEP':
        line = file.readline().strip()
        time = int(line)
      elif item == 'NUMBER':
        line = file.readline().strip()
        nAtoms = int(line)
      elif item == 'BOX':
        box = np.fromfile(file, float, 6, ' ')
        box = box.reshape((3,2))
      elif item == 'ATOMS':
        atoms = np.fromfile(file, float, nAtoms*nCols, ' ')
        atoms = atoms.reshape((nAtoms, nCols))
        break
      else:
        continue
    line = file.readline().strip()

  return time, box, atoms

def readFrame(file,nCols=4):
  line = file.readline().strip()
  nAtoms = int(line)

  line = file.readline()
  time = int(line.split()[-1])

  atoms = np.fromfile(file, float, nAtoms*nCols, ' ')
  atoms = atoms.reshape((nAtoms,nCols))

  return time, atoms

def write_xyz(filename, atoms, time=0, mode='w'):
  with open(filename+'.xyz',mode) as otp:
    otp.write('%d\nAtoms. Timestep: %d\n' % (len(atoms),time))
    try:
      [otp.write('%s %s %s %s\n' % tuple(line)) for line in atoms]
    except TypeError:
      [otp.write('%d %f %f %f\n' % tuple(line)) for line in atoms]

def write_conf(filename,
               atoms,
               bonds=[],
               title='Entangled polymer simulation with chains of N = 250\n',
               types=[1,0], # [atomtypes,bondtypes]
               box=[[-50,50],[-50,50],[0,100]],
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
