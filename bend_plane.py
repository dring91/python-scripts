from __future__ import division
from sys import exit, argv
from conf_tools import *
import numpy as np
import argparse

MODE = 'r'

def main():
  
  parser = argparse.ArgumentParser()
  parser.add_argument('-i','--input')
  parser.add_argument('-o','--output')
  parser.add_argument('-f','--formats',nargs='+',choices={'lammps','conf','xyz'},
                      default='conf')
  parser.add_argument('--line',nargs=4,help='axis normal loc angle')
  args = parser.parse_args()

  axis, normal, loc, angle = args.line
  loc, angle = float(loc), float(angle) * np.pi / 180

  with open(args.input, 'r') as file:
    box, atoms, _ = readConf(file)

  axes = {'x':0,'y':1,'z':2}
  const_dims = {axis, normal}
  divider = (axes.viewkeys() ^ const_dims).pop()
  surface = atoms[:,3:6]
  parallel = surface[:,axes[divider]]
  perpendicular = surface[:,axes[normal]]
  parallel += (parallel * (parallel < loc) - loc) * np.cos(angle)
  perpendicular += perpendicular * (perpendicular < loc) * np.sin(angle)
  surface[:,axes[divider]], surface[:,axes[normal]] = parallel, perpendicular
  atoms[:,3:6] = surface
  header = 'bent plane'
  print(np.cos(angle),np.sin(angle))

  # generate input files for lammps and VMD
  if 'conf' in args.formats:
    write_conf(args.output+'_out', atoms, title=header, box=box)
  if 'xyz' in args.formats:
    write_xyz(args.output+'_out', atoms[:,2:6])
  if 'lammps' in args.formats:
    write_traj(args.output+'_out', np.delete(atoms[:,:6],1,axis=1).astype(float), box, mode='w')

if __name__ == '__main__':
  main()
