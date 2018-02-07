import numpy as np

def readConf(file, cast=True):

  box = []
  atoms = []
  velocities = []
  bonds = []
  angles = []
  header = 'header'
  types = {}
  for line in file:
    L = line.split()
    if len(L) > 0 and ' '.join(L[1:3]) == 'atom types':
      types['atoms'] = int(L[0])
    if len(L) > 0 and ' '.join(L[1:3]) == 'bond types':
      types['bonds'] = int(L[0])
    if len(L) > 0 and ' '.join(L[1:3]) == 'angle types':
      types['angles'] = int(L[0])
    if len(L) > 0 and L[0] in ['Atoms','Velocities','Bonds','Angles']:
      header = L[0]
    if len(L) > 0 and L[-1] in ['xhi','yhi','zhi']:
      box.append(L[:2])
    elif len(L) > 4 and header == 'Atoms':
      atoms.append(L)
    elif len(L) > 4 and header == 'Velocities':
      velocities.append(L)
    elif len(L) > 3 and header == 'Bonds':
      bonds.append(L)
    elif len(L) > 3 and header == 'Angles':
      angles.append(L)

  if cast:
    box = np.array(box,dtype=float)
    atoms = np.array(atoms,dtype=float)
    velocities = np.array(velocities,dtype=float)
    bonds = np.array(bonds,dtype=int)
    angles = np.array(angles,dtype=int)
      
  return types, box, atoms, velocities, bonds, angles

def read_sections(file, sections):

  ## initialize a dictionary of all the types in the section
  types = OrderedDict()
  box = []
  config = OrderedDict((section,[]) for section in sections)
  ## add header comments to config and remember to 'pop' them before writing sections
  config['comments'] = []
  header,comment = 'header','comment'
  for line in file:
    L = line.split()
    if len(L) > 2 and L[2] == 'types':
      types[L[1]] = L[0]
    if len(L) > 3 and L[3] in ['xhi','yhi','zhi']:
      box.append(L[:2])
    if len(L) > 0 and L[0] in sections:
      header,comment = L[0],L[1:]
      config['comments'].append(comment) 
    elif len(L) > 0 and ' '.join(L[:2]) in sections:
      header,comment = ' '.join(L[:2]),L[2:]
      config['comments'].append(comment) 
    elif len(L) > 1 and header in config:
      config[header].append(L)

  return types, box, config

def tupleConf(file, cast=True):

  atoms = []
  box = []
  bonds = []
  header = 'header'
  for line in file:
    L = line.split()
    if len(L) > 0 and L[0] in ['Atoms','Bonds']:
      header = L[0]
    if len(L) > 0 and L[-1] in ['xhi','yhi','zhi']:
      box.append(tuple(L[:2]))
    elif len(L) > 2 and header == 'Atoms':
      atoms.append(tuple(L))
    elif len(L) > 2 and header == 'Bonds':
      bonds.append(tuple(L))

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
        nCols = len(L[2:])
        atoms = np.fromfile(file, float, nAtoms*nCols, ' ')
        atoms = atoms.reshape((nAtoms, nCols))
        break
      else:
        continue
    line = file.readline().strip()

  return time, box, atoms

def iTrj(file, nCols=5):

  line = file.readline()
  while (line != ''):
    # read timestep
    line = file.readline().strip()
    time = int(line)
    # read number of atoms
    line = file.readline()
    line = file.readline().strip()
    nAtoms = int(line)
    # read box size
    line = file.readline()
    box = np.fromfile(file, float, 6, ' ')
    box = box.reshape((3,2))
    # read atoms
    line = file.readline()
    atoms = np.fromfile(file, float, nAtoms*nCols, ' ')
    atoms = atoms.reshape((nAtoms, nCols))
    line = file.readline()

    yield time, nAtoms, box, atoms

def readFrame(file,nCols=4):
  line = file.readline().strip()
  nAtoms = int(line)

  line = file.readline()
  time = int(line.split()[-1])

  atoms = np.fromfile(file, float, nAtoms*nCols, ' ')
  atoms = atoms.reshape((nAtoms,nCols))

  return time, atoms

def write_traj(filename, atoms, box, time=0, mode='w', scaled=False, bnds=('pp','pp','pp')):
  if scaled: 
    atoms[:,2:5] = (atoms[:,2:5] - box[:,0]) / (box[:,1] - box[:,0])
    x,y,z = 'xs','ys','zs'
  else:
    x,y,z = 'x', 'y', 'z'
  atoms = np.nan_to_num(atoms)
  with open(filename+'.lammpstrj',mode) as file:
    file.write("ITEM: TIMESTEP\n")
    file.write(str(time)+'\n')
    file.write("ITEM: NUMBER OF ATOMS\n")
    file.write(str(len(atoms))+'\n')
    file.write("ITEM: BOX BOUNDS {} {} {}\n".format(*bnds))
  with open(filename+'.lammpstrj','ab') as file:
    np.savetxt(file, box, fmt="%.6f")
  with open(filename+'.lammpstrj','a') as file:
    file.write("ITEM: ATOMS id type {} {} {}\n".format(x,y,z))
  with open(filename+'.lammpstrj','ab') as file:
    np.savetxt(file, atoms, fmt="%d %d %f %f %f")

def write_xyz(filename, atoms, time=0, mode='w'):
  with open(filename+'.xyz',mode) as otp:
    otp.write('%d\nAtoms. Timestep: %d\n' % (len(atoms),time))
    try:
      [otp.write('%s %s %s %s\n' % tuple(line)) for line in atoms]
    except TypeError:
      [otp.write('%d %f %f %f\n' % tuple(line)) for line in atoms]

def write_sections(filename,
                   types=OrderedDict([("Atoms",1)]),
                   box=[['-1','1'],['-1','1'],['-1','1']],
                   config=OrderedDict([('Masses',[1]),('Atoms',[1,1,1,0.0,0.0,0.0,0,0,0])]),
                   title=''):
  
  with open(filename+'.conf','w') as file:  
    # write title line and skip a line
    file.write(title+'\n\n')
    
    # write number of atoms and number of bonds and skip a line
    file.write(str(len(config['Atoms']))+' atoms\n')
    if "bond" in types: file.write(str(len(config['Bonds']))+' bonds\n')
    if "angle" in types: file.write(str(len(config['Angles']))+' angles\n')
    file.write('\n')
    
    # write types and skip a line
    for name, num in types.items():
      file.write(' '.join([num, name, 'types\n']))
    file.write('\n')
    
    # write box dimensions and skip a line
    dims = 'xyz'
    [file.write('{0} {1} {2}lo {2}hi\n'.format(lo,hi,dim)) for dim,(lo,hi) in zip(dims,box)]
    file.write('\n')

    # pop off header comments
    comments = config.pop('comments')

    # construct body with sections
    for comment,(section,data) in zip(comments,config.items()):
      print(section)
      file.write('{} {}\n'.format(section,' '.join(comment)))
      file.write('\n')
      [file.write(" ".join(line)+'\n') for line in data]
      file.write('\n')

def write_conf(filename,
               atoms,
               bonds=[],
               angles=[],
               box=[[-50,50],[-50,50],[0,100]],
               types={"atoms":1, "bonds":0, "angles":0},
               masses=[1], 
               title=''):
  
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
    if types["bonds"] > 0: otp.write(str(len(bonds))+' bonds\n')
    if types["angles"] > 0: otp.write(str(len(angles))+' angles\n')

    # skip line
    otp.write('\n')
    
    # write types
    otp.write(str(types["atoms"])+' atom types\n')
    if types["bonds"] > 0: otp.write(str(types["bonds"])+' bond types\n')
    if types["angles"] > 0: otp.write(str(types["angles"])+' angle types\n')
    
    # skip line
    otp.write('\n')
    
    # write box dimensions
    dims = 'xyz'
    [otp.write('{1:.4f} {2:.4f} {0}lo {0}hi\n'.format(*dim)) for dim in zip(dims,*box.T)]
    
    # skip line
    otp.write('\n')

    # write masses
    otp.write('Masses\n')
    otp.write('\n')
    [otp.write('{} {}\n'.format(i+1,mass)) for i,mass in enumerate(masses)]

    # skip line    
    otp.write('\n')
    
    # write atoms
    otp.write('Atoms\n')
    otp.write('\n')
  with open(filename+'.conf','ab') as otp: np.savetxt(otp, atoms, fmt='%d %d %d %f %f %f %d %d %d')
    #[otp.write('{} {} {} {} {} {}\n'.format(*line)) for line in atoms]

  if types["bonds"] > 0:
    with open(filename+'.conf','a') as otp:  
      # skip line    
      otp.write('\n')
      
      # write bonds
      otp.write('Bonds\n')
      otp.write('\n')
    with open(filename+'.conf','ab') as otp:  
      np.savetxt(otp, bonds, fmt='%d')
      #[otp.write('{} {} {} {}\n'.format(*line)) for line in bonds]
  if types["angles"] > 0:
    with open(filename+'.conf','a') as otp:  
      # skip line    
      otp.write('\n')
      
      # write bonds
      otp.write('Angles\n')
      otp.write('\n')
    with open(filename+'.conf','ab') as otp:  
      np.savetxt(otp, angles, fmt='%d')
      #[otp.write('{} {} {} {}\n'.format(*line)) for line in bonds]
