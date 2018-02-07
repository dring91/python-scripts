from sys import argv, exit
import argparse
from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np
from copy import copy
from random import choice

def read_sections(file, sections):

  ## initialize a dictionary of all the types in the section
  types = OrderedDict()
  box = []
  config = OrderedDict((section,[]) for section in sections)
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

def create_angles(chains,mon,offset):
  angles = [[m*(mon-2)+n+1,1,m*mon+n+offset,m*mon+n+1+offset,m*mon+n+2+offset] for m in range(chains) for n in range(mon-2)]
  angles = [[str(angle) for angle in row] for row in angles]
  return angles

def main():

  """ Add, modify, or subtract section to configuration file """

  ## List of the possible sections that can be in a config file
  sections = ['Masses','Pair Coeffs','Bond Coeffs','Angle Coeffs','Atoms','Velocities','Bonds','Angles']
  sections = ['Masses','Atoms','Bonds']

  # parse input
  parser = argparse.ArgumentParser()
  parser.add_argument("-i","--input",help="name of output file sans extension")
  parser.add_argument("-o","--output",help="name of output file sans extension")
  parser.add_argument("-c","--chains",type=int,help="Number of chains in melt")
  parser.add_argument("-m","--mon",type=int,help="Number of monomers to a chain")
  parser.add_argument("--sort",action='store_true')
  parser.add_argument("--add",nargs='+',help="sections to add to the config file")
  parser.add_argument("--sub",nargs='+',help="sections to subtract from the config file")
  parser.add_argument("--mod",nargs='+',help="sections to modify in the config file")
  parser.add_argument("--title",required=True)
  args = parser.parse_args()

  with open(args.input, 'r') as file:
    types, box, config = read_sections(file, sections)

  if args.sort is not None:
    config['Atoms'].sort(key=lambda x:int(x[0]))
    config['Bonds'].sort(key=lambda x:int(x[0]))
    #config['Velocities'].sort(key=lambda x:int(x[0]))
    #config['Angles'].sort(key=lambda x:int(x[0]))

  if args.add is not None:
    if 'Angles' in args.add:
      config['Bonds'] = [[str(i+1)]+bond[1:] for i,bond in enumerate(sorted(config['Bonds'],key=lambda x:int(x[2])))]
      #config['Velocities'] = [atom for atom in sorted(config['Velocities'],key=lambda x:int(x[0]))]
      offset = min(int(atom[0]) for atom in config['Atoms'] if atom[2] == '2')
      angles = [(args.add[0],create_angles(args.chains,args.mon,offset))]
      config.update(angles)
      types.update([(args.add[0].lower()[:-1],'1')])
      config['comments'].append([])

  if args.mod is not None:
    if 'unwrap' in args.mod:
      lengths = [float(hi) - float(lo) for lo,hi in box] 
      config['Atoms'] = [atom for atom in config['Atoms']]
      config['Atoms'] = [atom[:3] + [str(float(x) + l*int(f)) for x,l,f in zip(atom[3:6],lengths,atom[6:9])] + ['0']*3 for atom in config['Atoms']]
    if 'numbering' in args.mod:
      ## adjust molecule numbering
      config['Atoms'] = [atom for atom in config['Atoms']]
      polymers = [atom for atom in config['Atoms'] if atom[2] == '2']
      ## change the molecule ids for bond/swap
      mid_chain = args.mon//2 + args.mon%2
      molecule = [m+1 for m in range(mid_chain)] + \
                 [m+1 for m in reversed(range(args.mon//2))]
      molecules = [str(m+1) for _ in range(args.chains) for m in molecule]
      mol = (str(2+mid_chain),str(3+mid_chain))
      atoms = [atom[:1] + ['1'] + atom[2:] for atom in config['Atoms'] if atom[2] == '1'] + \
              [atom[:1] + [mol] + atom[2:] for atom,mol in zip(polymers,molecules)] + \
              [atom[:1] + [mol[0]] + atom[2:] for atom in config['Atoms'] if atom[2] == '3'] + \
              [atom[:1] + [mol[1]] + atom[2:] for atom in config['Atoms'] if atom[2] == '4']
      ## renumber atom ids
      ids = [str(i+1) for i,atom in enumerate(atoms)]
      diffs = {old[0]:new for old,new in zip(atoms,ids)}
      atoms = [[str(i+1)] + atom[1:] for i,atom in enumerate(atoms)]
      #config['Velocities'] = [[id] + vel[1:] for id,vel in zip(ids,config['Velocities'])]
      config['Bonds'] = [bond[:2] + [diffs[bond[2]],diffs[bond[3]]] for bond in config['Bonds']]
      config['Angles'] = [bond[:2] + [diffs[bond[2]],diffs[bond[3]],diffs[bond[4]]] 
                                                       for bond in config['Angles']]
      bonds = {atom for bond in config['Bonds'] for atom in bond[2:]}
      ## use set subtraction to check that the bond atoms are consistent
      ids = set(ids)
      print(len(ids-bonds))
      config['Atoms'] = atoms
    if 'atom' in args.mod:
      ## select atoms with type=4
      plate = [atom for atom in config['Atoms'] if atom[2] == '3']
      cyl = [atom for atom in config['Atoms'] if atom[2] == '4']
      atoms = [atom for atom in config['Atoms'] if atom[2] != '4' and atom[2] != '3']
      ## iterate through atoms and randomly assign new types to each atom
      platetypes, cyltypes = ['3','4','5','6','7'], ['8','9','10','11','12']
      plate = [row[:2] + [choice(platetypes)] + row[3:] for row in plate]
      cyl = [row[:2] + [choice(cyltypes)] + row[3:] for row in cyl]
      ## update config
      atoms = [('Atoms',plate + cyl + atoms)]
      config.update(atoms)
      ## update atom types
      types.update([('atom','12')])
      ## update masses
      masses = [[str(i+1),'1'] for i in range(int(types['atom']))]
      config.update([('Masses',masses)])
    if 'bond' in args.mod:
      ## select only polymers
      polymers = [atom for atom in config['Atoms'] if atom[2] == '2']
      polymers.sort(key=lambda x:int(x[0]))
      atoms = [atom for atom in config['Atoms'] if atom[2] != '2']

      ## change the molecule ids for bond/swap
      molecule = [m+1 for m in range(args.mon//2+args.mon%2)] + \
                 [m+1 for m in reversed(range(args.mon//2))]
      molecules = [str(m) for _ in range(args.chains) for m in molecule]
      polymers = [row[:1] + [mol] + row[2:] for row,mol in zip(polymers,molecules)]
      config['Atoms'] = polymers + atoms

      ## remove bonds connected to the 25th and 26th atoms
      chain_ends = [atom[0] for atom in polymers if atom[1] == '1']
      bonds = [(a,b) for _,_,a,b in sorted(config['Bonds'],key=lambda x:int(x[0]))]
      extra_bonds = list(zip(chain_ends[1:-2:4],chain_ends[2:-1:4]))
      indices = [bonds.index(pair) for pair in extra_bonds]
      for index in sorted(indices,reverse=True): config['Bonds'].pop(index)

      ## remove angles with atoms 24,25,26 and 25,26,27
      chain_ends = [atom[0] for atom in polymers if atom[1] == '1' or atom[1] == '2']
      angles = [(a,b,c) for _,_,a,b,c in sorted(config['Angles'],key=lambda x:int(x[0]))]
      extra_angles = list(zip(chain_ends[2::8],chain_ends[3::8],chain_ends[4::8])) + \
                     list(zip(chain_ends[3::8],chain_ends[4::8],chain_ends[5::8]))
      indices = [angles.index(triple) for triple in extra_angles]
      for index in sorted(indices,reverse=True): config['Angles'].pop(index)

  if args.sub is not None:
    for section in args.sub:
      if section in config:
        del config[section]

  # output coordinates
  write_conf(args.output, types, box, config, title=args.title)
  
if __name__ == '__main__':
  main()
