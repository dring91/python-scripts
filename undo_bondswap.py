import numpy as np
from sys import argv, exit
from copy import copy

from conf_tools import readConf, write_conf, write_traj
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from numpy.linalg import norm

def check_bond_lengths(coords, bonds, plot=True):
  ## calculate the bond vector
  diffs = coords[bonds[:,2]-1] - coords[bonds[:,3]-1]
  ## calculate bond length
  dists = np.sqrt(diffs[:,0]**2+diffs[:,1]**2+diffs[:,2]**2)
  ## generate histograms of bond length and vector components
  fx, x = np.histogram(np.abs(diffs[:,0]), bins='auto')
  fy, y = np.histogram(np.abs(diffs[:,1]), bins='auto')
  fz, z = np.histogram(np.abs(diffs[:,2]), bins='auto')
  fr, r = np.histogram(dists, bins='auto') # should be symmetric and centered at 1
  
  ## plot histograms of bonds
  if plot:
    fig, (comps, total) = plt.subplots(ncols=2)
    comps.plot(x[:-1]+np.diff(x)/2,fx,label='x') 
    comps.plot(y[:-1]+np.diff(y)/2,fy,label='y') 
    comps.plot(z[:-1]+np.diff(z)/2,fz,label='z') 
    total.semilogy(r[:-1]+np.diff(r)/2,fr,label='r') 
    comps.legend()
    total.legend()
    fig.tight_layout()
    plt.show()

def rearrange_molecules(atoms, bonds):

  nChains = 1000
  nMon = 100
  new_atoms = np.zeros((nChains, nMon, 9))
  for n in range(nMon-1):
    heads = bonds[np.isin(bonds[:,2],bonds[:,3],invert=True)]
    bonds = bonds[np.isin(bonds[:,2],bonds[:,3],invert=False)]
    new_atoms[:,n] = atoms[heads[:,2]-1]

  return new_atoms
 
def main():
  
  parser = ArgumentParser()
  parser.add_argument('-i','--input')
  args = parser.parse_args()

  with open(args.input, 'r') as file:
    types, box, atoms, bonds = readConf(file)

  ## sort atoms by atomid
  atoms = atoms[atoms[:,0].argsort()]
  ## unwrap atoms using image flags
  atoms[:,3:6] = atoms[:,3:6] + atoms[:,6:9] * (box[:,1] - box[:,0])
  ## calculate bond lengths to confirm proper unwrapping and bonding of molecules
  check_bond_lengths(atoms[:,3:6], bonds, plot=False)
  ## rearrange bonds
  atoms = rearrange_molecules(atoms, bonds)
  ## Calculate end-to-end distances
  diffs = np.diff(atoms[:,:,3:6],axis=1)
  dists = diffs.sum(1)
  print(norm(diffs,axis=1), norm(dists,axis=1))
    
if __name__ == '__main__':
  main()
