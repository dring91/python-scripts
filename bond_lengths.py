import numpy as np
from conf_tools import *
from sys import argv
from math import sqrt

def main():

  try: 
    filename = argv[1]
  except IndexError:
    print 'Script requires a filename'

  types = ['1']
  with open(filename+'.conf','r') as file:
    box, atoms, bonds = readConf(file, types)

  nMon = 50
  limits = np.array(box, dtype=float)
  atoms = np.array(atoms)

  # flags = atoms[:,6:].astype(int)
  chains = atoms[:,3:6].astype(float)
  # chains = chains + (limits[:,1]-limits[:,0])*flags
  lengths = [[i+1,sqrt(sum([(chains[i,j]-chains[i+1,j])**2 for j in range(3)]))] 
                                                      for i in range(len(chains)-1)]

  with open(filename+'_lengths','w') as otp:
    [otp.write("%d %f\n" % tuple(line)) for line in lengths if line[0] % nMon]

if __name__ == '__main__':
  main()
