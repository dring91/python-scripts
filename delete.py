import numpy as np
import argparse
from conf_tools import *

def main():
  '''deletes/selects atoms/bonds from a configuration and writes them to a new file'''

  parser = argparse.ArgumentParser()
  parser.add_argument('-i','--input')
  parser.add_argument('-o','--output')
  parser.add_argument('-t','--types',nargs='+',default='')
  parser.add_argument('-b','--bonds',nargs='+',default='')
  parser.add_argument('-d','--description',default='')
  args = parser.parse_args()
  
  with open(args.input) as file:
    box, atoms, bonds = readConf(file, args.types)

  atoms = atoms[np.argsort(atoms[:,0])]

  n = len(args.types)
  m = len(args.bonds)
  write_conf(args.output, atoms.astype('|S16'), bonds, box, {"atoms":n, "bonds":m}, [1]*n, args.description)
  write_traj(args.output, np.delete(atoms[:,:6],1,1), box)

if __name__ == '__main__':
  main()
