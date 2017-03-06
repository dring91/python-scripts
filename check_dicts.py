import numpy as np
import argparse

from conf_tools import readConf

def main():
  
  # get commandline arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("-i","--input",help='input file to renumber')
  args = parser.parse_args()
  
  # read input file
  with open(args.input,'r') as file:
    box, frame, bonds = readConf(file)

  ## Dictionaries to create:
  # 1. Dictionary of atoms
  atoms = {atom[0]:atom[3:] for atom in frame} # np.delete(frame,[1,2],1)
  # 2. Dictionary of molecules
  molIds = set(frame[:,1])
  molecules = {}
  for mol in molIds:
    molecules[mol] = frame[:,0][frame[:,1] == mol]
  # 3. Dictionary of types
  types = set(frame[:,2])
  atomtypes = {}
  for type in types:
    atomtypes[type] = frame[:,:1][frame[:,2] == type]
  # 4. Dictionary of bonds
  ends = bonds[:,2][np.in1d(bonds[:,2],bonds[:,3],invert=True)]
  bonds = {bond[2]:bond[3] for bond in bonds}

  ## Display dictionaries:
  print {k:v for k,v in atoms.iteritems() if k > 5 and k < 10}

  ## Tests to perform:
  # Are all the molecule bonded correctly?
  #   1. Are the chains captured properly in the dictionary?
  #   2. Is there agreement between molecules and bonds?

if __name__ == "__main__":
  main()
