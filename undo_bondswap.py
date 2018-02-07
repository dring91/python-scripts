import numpy as np
from sys import argv, exit
from copy import copy
from argparse import ArgumentParser
from collections import OrderedDict

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

def getEnds(atoms, bonds):
  # find them by counting up the number of times an atom occurs in bonds
  atom_count = {key:0 for key,value in atoms.items()}
  bond_atoms = [item for key,value in bonds.items() 
                     for item in (key, value)]
  for atom in bond_atoms:
    atom_count[atom] += 1
  ends = [key for (key,value) in atom_count.items() if value == 1]

  return ends

def unshuffleAtoms(atoms, bonds, ends, nMon):
  # loop over bond ends array
  # build a new list of atoms by picking a chain end and adding atoms
    # when an atom is added, remove it from the atoms list
  # if chainLength is reached and the last atom is an end, you're good
  # repeat until all bond ends are used up

  unshuffled = []
  # unshuffle all the 'forwards' chains and then flip the 'backwards' ones and unshuffle the remaining
  for tail in ends:
    if tail in bonds.keys():
      atom = tail
      unshuffled.append([atom] + atoms[atom]) 
      for n in range(nMon-1):
        atom = bonds.pop(atom)
        unshuffled.append([atom] + atoms.pop(atom)) 
    elif tail in bonds.values():
      # print 'Value in dict reversed'
      continue

  return unshuffled

def write_conf(filename,
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

def main():
  
  ## List of the possible sections that can be in a config file
  sections = ['Masses','Pair Coeffs','Bond Coeffs','Angle Coeffs','Atoms','Velocities','Bonds','Angles']
  sections = ['Masses','Atoms','Bonds','Angles']

  # parse commandline args
  parser = ArgumentParser()
  parser.add_argument('-i','--input')
  parser.add_argument('-l','--chainlength',type=int)
  parser.add_argument('--polymer')
  args = parser.parse_args()

  with open(args.input+'.conf', 'r') as file: types, box, config = read_sections(file, sections)

  # unwrap atom coordinates using provided image flags
  side_lengths = [float(row[1]) - float(row[0]) for row in box]
  coordinates = [list(map(float,row[3:6])) for row in config['Atoms']]
  image_flags = [list(map(int,row[6:])) for row in config['Atoms']]
  coordinates = [[str(x+l*f) for x,l,f in zip(particle,side_lengths,image)] 
                             for particle,image in zip(coordinates,image_flags)]
  image_flags = [['0','0','0'] for row in image_flags]
  config['Atoms'] = [ids[:3] + coords + flags for ids,coords,flags in zip(config['Atoms'],coordinates,image_flags)]

  # search for all the chain ends in bonds and store in new array
  atoms = {row[0]:row[1:] for row in config['Atoms'] if row[2] == args.polymer}
  bonds = {row[2]:row[3] for row in config['Bonds']}
  ends = getEnds(atoms, bonds)

  # unshuffle atoms
  unshuffled = unshuffleAtoms(copy(atoms), copy(bonds), ends, args.chainlength)
  atoms = unshuffled + [row for row in config['Atoms'] if row[2] != args.polymer]
  config['Atoms'] = atoms

  # write the new atoms list and the bonds to a file
  write_conf(args.input+'_unshuffled', types, box, config, title='unshuffled configuration')

if __name__ == '__main__':
  main()
