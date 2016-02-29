import numpy as np
from sys import argv, exit
from getopt import *
import matplotlib.pyplot as plt
from conf_tools import *

def writeHeader(filename,header):
  # Implicitly assumes a 1 line header
  with open(filename, 'w') as file:
    file.write(header)   
  
def writeChains(filename,chains,time,mode):
  with open(filename, mode) as file:
    [file.write('%d %f %f %f %f\n' % (time, com[0], com[1], com[2], com[3])) for com in chains]
    file.write('\n\n')

def getArgs(argv):

  try:
    opts, args = getopt(argv[1:],'i:o:n:')
  except GetoptError:
    print 'randomWalk.py -i <xyz file> -o <filename> -n <nFrames>'
    exit(2)
  for opt, arg in opts:
    if opt == '-i':
      inFile = arg
    elif opt == '-o':
      outFile = arg
    elif opt == '-n':
      nSteps = arg
 
  return inFile, outFile, int(nSteps)
  
def main():

  trajFile, outFile, nFrames = getArgs(argv)

  title = '# polymer chain trajectories\n'+'# time  COM_z   RG\n'
  writeHeader(outFile,title)

  mode = 'a'

  with open(trajFile, 'r') as file:
    for n in range(nFrames):
      # read data
      t, frame = readFrame(file)

      # filter and reshape data
      polymers = frame[frame[:,0] == 1][:,1:]
      polymers = polymers.reshape((10,-1,3),order='F')

      # Calculate COM and RG
      COM = polymers.mean(0)
      RG = (polymers - COM)**2
      RG = RG.sum(2)
      RG = np.sqrt(RG.mean(0))

      # output
      if len(RG) != 0:
        # write COM, RG and surface atoms
        output = np.zeros((len(RG),4))
        output[:,:3] = COM
        output[:,3] = RG
        writeChains(outFile, output, t, mode)
        print('%d Chains written at timestep: %d' % (len(output),t))
      else:
        print 'No observed infiltration'
          
  print('Finished analyzing trajectory')

if __name__ == "__main__":  
  main()
