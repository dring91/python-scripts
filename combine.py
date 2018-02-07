from sys import argv, exit
import numpy as np
from argparse import ArgumentParser
from functools import reduce
from collections import OrderedDict

from pbc_tools import *
from conf_tools import readConf, write_conf, write_traj, write_xyz

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

#def write_conf(filename,
#               types=OrderedDict([("Atoms",1)]),
#               box=[['-1','1'],['-1','1'],['-1','1']],
#               config=OrderedDict([('Masses',[1]),('Atoms',[1,1,1,0.0,0.0,0.0,0,0,0])]),
#               title=''):
#  
#  with open(filename+'.conf','w') as file:  
#    # write title line and skip a line
#    file.write(title+'\n\n')
#    
#    # write number of atoms and number of bonds and skip a line
#    file.write(str(len(config['Atoms']))+' atoms\n')
#    if "bond" in types: file.write(str(len(config['Bonds']))+' bonds\n')
#    if "angle" in types: file.write(str(len(config['Angles']))+' angles\n')
#    file.write('\n')
#    
#    # write types and skip a line
#    for name, num in types.items():
#      file.write(' '.join([num, name, 'types\n']))
#    file.write('\n')
#    
#    # write box dimensions and skip a line
#    dims = 'xyz'
#    [file.write('{0} {1} {2}lo {2}hi\n'.format(lo,hi,dim)) for dim,(lo,hi) in zip(dims,box)]
#    file.write('\n')
#
#    # pop off header comments
#    comments = config.pop('comments')
#
#    # construct body with sections
#    for comment,(section,data) in zip(comments,config.items()):
#      print(section)
#      file.write('{} {}\n'.format(section,' '.join(comment)))
#      file.write('\n')
#      [file.write(" ".join(line)+'\n') for line in data]
#      file.write('\n')

def main():

  parser = ArgumentParser(description="Handle options for main program")
  parser.add_argument("-f","--file",action="append",nargs='+')
  parser.add_argument("-o", "--output", help="output file name")
  parser.add_argument("-d", "--descrip", help="description for output file header", 
                      default="")
  parser.add_argument("--formats", nargs='+', help="output formats",
                      choices=['xyz','lammps','conf'], default='conf')
  parser.add_argument("--no-image-flags",dest="entries",action="store_const",const=6,
                      default=9,help="switch to adjust number of entries for data flags")
  parser.add_argument("--box",action="append",nargs=2)
  args = parser.parse_args()

  subpar = ArgumentParser(description="parse input for each file")
  subpar.add_argument("file_name")
  subpar.add_argument("type",type=int)
  subpar.add_argument("shift",nargs=3,type=float)

  optpar = ArgumentParser(description="parse box options to expand dimensions")
  optpar.add_argument("wall",choices={"xlo","xhi","ylo","yhi","zlo","zhi"})
  optpar.add_argument("coord",type=float)

  # initialize list of layers
  boxes = []
  atomic = []
  bonded = []
  angled = []
  type_list = []
  atomtypes = set()
  bondtypes = 1
  atom_ids = 0
  mol_ids = 0
  # initialize parser dictionaries
  dims = {'x':0,'y':1,'z':2}
  faces = {'lo':0,'hi':1}
  # parse file information 
  for file_entry in args.file:
    file_data = subpar.parse_args(file_entry)
    # read data files
    with open(file_data.file_name, 'r') as file:
      types, box, atoms, velocities, bonds, angles = readConf(file)
      type_list.append(types)
      box += np.array(file_data.shift)[:,np.newaxis]
      boxes.append(box)
      atoms[:,0] += atom_ids
      atoms[:,1] += mol_ids
      atoms[:,2] = file_data.type
      atoms[:,3:6] += file_data.shift
      atomic.append(atoms[atoms[:,0].argsort()])
      atomtypes.add(file_data.type)
      if len(bonds) > 0:
        bonds = bonds[bonds[:,0].argsort()]
        bonds[:,2:] = bonds[:,2:] + atom_ids
        bonded.append(bonds)
      if len(angles) > 0:
        angles = angles[angles[:,0].argsort()]
        angles[:,2:] = angles[:,2:] + atom_ids
        angled.append(angles)
      atom_ids += atoms[:,0].max()
      mol_ids +=  atoms[:,1].max()

  def merge(d1, d2): return { **d1, **d2 }
  types = reduce(merge, type_list)
  if types.get('angles') is None: types.update(angles=0)

  boxes = np.stack(boxes,axis=2)
  box[:,0], box[:,1] = boxes[:,0,:].min(1), boxes[:,1,:].max(1)
  if args.box is not None:
    for box_data in args.box:
      opt = optpar.parse_args(box_data)
      dim, face = dims[opt.wall[0]], faces[opt.wall[1:]]
      box[dim,face] = opt.coord
    
  atoms = np.concatenate(atomic,axis=0)
  bonds = np.concatenate(bonded,axis=0)
  bonds = bonds[bonds[:,2].argsort()]
  bonds[:,0] = np.arange(len(bonds))+1
  if len(angled) > 0:
    angles = np.concatenate(angled,axis=0)
    angles = angles[angles[:,2].argsort()]
    angles[:,0] = np.arange(len(angles))+1

  # write output files
  types["atoms"] = len(atomtypes)
  if 'xyz' in args.formats:
    write_xyz(args.output, [line[2:6] for line in atoms])
  if 'conf' in args.formats:
    write_conf(args.output, atoms, bonds, angles, box, types, [1] * types["atoms"], args.descrip)
  if 'lammps' in args.formats:
    atoms[:,3:6] += atoms[:,6:9] * (box[:,1] - box[:,0])
    write_traj(args.output, np.delete(atoms[:,:6],1,1), box, mode='w')

if __name__ == "__main__":
  main()
