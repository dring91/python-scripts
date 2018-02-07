from sys import argv, exit
from argparse import ArgumentParser
from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np
from copy import copy
from random import choice
from math import sqrt
from itertools import islice,chain

def read_sections(file, sections):

  ## initialize a dictionary of all the types in the section
  types, box, config = OrderedDict(), [], OrderedDict((section,[]) for section in sections)
  ## add header comments to config and remember to 'pop' them before writing sections
  config['comments'] = []
  header, comment = 'header', 'comment'
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

def bond_lengths(config):
  pass

def main():

  """ check that the conf file doesn't contain inconsistencies """

  ## List of the possible sections that can be in a config file
  sections = ['Masses','Pair Coeffs','Bond Coeffs','Atoms','Velocities','Bonds','Angles']

  # parse input
  parser = ArgumentParser()
  parser.add_argument("-i","--input",help="name of output file sans extension")
  parser.add_argument("-o","--output",help="name of output file sans extension")
  parser.add_argument("-c","--chains",type=int,help="Number of chains in melt")
  parser.add_argument("-m","--mon",type=int,help="Number of monomers to a chain")
  args = parser.parse_args()

  with open(args.input, 'r') as file:
    types, box, config = read_sections(file, sections)

  """ What quantities need to be checked?
      1. number of atoms, bonds, masses, etc. should all be consistent
      2. atom, bond, and angle ids should be unique
      3. mol should be unique and numbered properly
      4. types should all agree
      5. chain statistics and bond topology should be consistent
      6. image flags should be consistent
      7. formatting should be good
  """

  ## unwrap atoms
  lengths = [float(hi) - float(lo) for lo, hi in box]
  config['Atoms'] = [[float(val) for val in row] for row in config['Atoms']]
  config['Atoms'].sort(key=lambda x:x[0])
  config['Atoms'] = [row[:3] + [u+l*f for u,l,f in zip(row[3:6],lengths,row[6:9])] 
                                                        for row in config['Atoms']]
  ## calculate bond lengths
  config['Bonds'] = [[int(val) for val in row] for row in config['Bonds']]
  bonds = [[config['Atoms'][bond[2]-1][3:],config['Atoms'][bond[3]-1][3:]] for bond in config['Bonds']]

  bond_lengths = [sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2) for (x,y,z),(x0,y0,z0) in bonds]
  mean = sum(bond_lengths)/len(bond_lengths)
  print(mean, sqrt(sum((length-mean)**2 for length in bond_lengths)/len(bond_lengths)))
  print("b_min: {} b_max: {}".format(min(bond_lengths),max(bond_lengths)))

  ## calculate chain end-to-end distances
  vecs = [[x0-x,y0-y,z0-z] for (x,y,z),(x0,y0,z0) in bonds]  
  bond_vectors = np.array(vecs)
  bond_vectors = bond_vectors.reshape((-1,args.mon-1,3))
  end_to_end = bond_vectors.sum(axis=1)
  end_to_end = np.linalg.norm(end_to_end,axis=1)
  print(end_to_end.mean())

  ## print out sets of particle types, numbers etc.
  ids = {atom[0] for atom in config['Atoms']}
  mol = {atom[1] for atom in config['Atoms']}
  typ = {atom[2] for atom in config['Atoms']}
  print(len(ids),min(ids),max(ids))
  print(len(mol),min(mol),max(mol))
  print(len(typ),min(typ),max(typ))

  ## create a dictionary of atom types and molecules
  mol_type_mapping = {atom[2]:set([]) for atom in config['Atoms']}
  for atom in config['Atoms']:
    mol_type_mapping[atom[2]].add(atom[1])

  for typ,mols in mol_type_mapping.items():
    print(typ, mols)

if __name__ == '__main__':
  main()
