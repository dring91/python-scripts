import numpy as np
from sys import argv, exit
from getopt import *
import matplotlib.pyplot as plt
from conf_tools import *

def writeHeader(filename,header):
  # Implicitly assumes a 1 line header
  with open(filename, 'w') as file:
    file.write(header)   
  
def writeChains(filename,chains,mode):
  with open(filename, mode) as file:
    [file.write('%f %f %f %f %f\n' % tuple(com)) for com in chains]
    file.write('\n\n')

def getArgs(argv):

  try:
    opts, args = getopt(argv[1:],'i:o:n:')
  except GetoptError:
    print 'randomWalk.py -i <xyz file> -o <filename>'
    exit(2)
  for opt, arg in opts:
    if opt == '-i':
      inFile = arg
    elif opt == '-o':
      outFile = arg
 
  return inFile, outFile
  
def main():

  trajFile, outFile = getArgs(argv)

  title = '# polymer chain trajectories\n'+'# time  COM_z   RG\n'
  writeHeader(outFile,title)

  mode = 'a'
  t = 0
  nMon = 25

  with open(trajFile, 'r') as file:
    # read data
    box, atoms, bonds = readConf(file,'1')

    box = np.array(box,dtype=float)
    bonds = np.array(bonds,dtype=int)
    atoms = np.array(atoms)
    info = np.array(atoms[:,:3],dtype=int)
    polymers = np.array(atoms[:,3:6],dtype=float)
    flags = np.array(atoms[:,6:],dtype=int)

    # reshape data
    polymers = polymers.reshape((nMon,-1,3),order='F')

    # Calculate COM and RG
    COM = polymers.mean(0)
    RG = (polymers - COM)**2
    RG = RG.sum(2)
    RG = np.sqrt(RG.mean(0))
    RE = polymers[-1] - polymers[0]
    RE = np.sqrt(RE[:,0]**2 + RE[:,1]**2 + RE[:,2]**2)

    # output
    if len(RG) != 0:
      # write COM, RG and surface atoms
      output = np.zeros((len(RG),5))
      output[:,:3] = COM
      output[:,3] = RG
      output[:,4] = RE
      writeChains(outFile, output, mode)
      # print('%d Chains written at timestep: %d' % (len(output),t))
    else:
      print 'No observed infiltration'
          
  print('Finished analyzing trajectory')

if __name__ == "__main__":  
  main()
