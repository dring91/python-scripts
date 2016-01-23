#!/usr/bin/python

from sys import argv, exit
from getopt import *
from writeAtoms import *
from nestedDicts import *
import numpy as np

def getArgs(argv):

  try:
    opts, args = getopt(argv[1:],'i:o:')
  except GetoptError:
    print 'randomWalk.py -i <filename> -o <filename>'
    exit(2)
  for opt, arg in opts:
    if opt == '-i':
      inFile = arg
    elif opt == '-o':
      outFile = arg
 
  return inFile, outFile
  
def readConf(file, atype):

  atoms = []
  box = []
  bonds = []
  header = 'header'
  for line in file:
    L = line.split()
    if len(L) > 0 and L[0] in set(['Atoms','Bonds']):
      header = L[0]
    if len(L) > 0 and L[-1] in set(['xhi','yhi','zhi']):
      box.append(L[:2])
    elif len(L) > 2 and L[2] in set(atype) and header == 'Atoms':
      atoms.append(L)
    elif len(L) > 2 and header == 'Bonds':
      bonds.append(L)
      
  return np.array(box,dtype=float), np.array(atoms,dtype=float), np.array(bonds,dtype=float)
  
def main():
  
  inFile, outFile = getArgs(argv)
  
  with open(inFile,'r') as file:
    box, atoms, bonds = readConf(file, ['3','2','1'])
  
  partMask = (atoms[:,2] == 3) 
  polyMask = (atoms[:,2] == 1)
  surfMask = (atoms[:,2] == 2)
  particles = (partMask*atoms.T).T
  polymers = (polyMask*atoms.T).T
  surface = (surfMask*atoms.T).T

  polyMin = polymers[:,0][polyMask].min()
  polymers[:,1] = (polymers[:,0] - polyMin + particles[:,1].max() + 1)*polyMask
  surface[:,1] = (polymers[:,1].max() + 1)*surfMask

  atoms = particles + polymers + surface
  
  title = 'unwrapped undersaturated configuration'
  types = [3,0]
  masses = [1,1,1]  
  write_conf(outFile,atoms,bonds,title,types,box,masses)
  write_xyz(outFile,atoms[:,2:])  
  
if __name__ == "__main__":
  main()
