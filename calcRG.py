import numpy as np
from sys import argv, exit
from getopt import *
import matplotlib.pyplot as plt
from conf_tools import *
from pbc_tools import *

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

def make2dHist(X,Y,limits):
  nBins = (50,50)
  hist = np.zeros((nBins[0]*nBins[1],4))
  binSize = ((limits[0,1] - limits[0,0]) / nBins[0], 
             (limits[1,1] - limits[1,0]) / nBins[1])
  xs = np.arange(nBins[0])
  ys = np.arange(nBins[1])
  hist[:,0] = np.repeat(xs, nBins[1])
  hist[:,1] = np.tile(ys, nBins[0])

  bins = ((X - limits[:,0])/binSize).astype(int)
  for i,bin in enumerate(hist[:,:2]):
    mask = (bins == bin)
    mask = mask[:,0]*mask[:,1]
    if mask.sum() != 0:
      hist[i,2] = Y[mask].mean()
      hist[i,3] = Y[mask].std()

  hist[:,0] = (hist[:,0] + 0.5) * binSize[0] + limits[0,0]
  hist[:,1] = (hist[:,1] + 0.5) * binSize[1] + limits[1,0]

  return hist
  
def main():

  trajFile, outFile, nFrames = getArgs(argv)

  title = '# polymer chain trajectories\n'+'# time  COM_z   RG\n'
  writeHeader(outFile,title)

  mode = 'a'
  nMon = 25

  with open(trajFile, 'r') as file:
    for n in range(nFrames):
      # read data
      t, frame = readFrame(file)
      # t, box, frame = readTrj(file)

      # filter and reshape data
      polymers = frame[frame[:,0] == 1][:,1:]
      # frame = frame[np.argsort(frame[:,0])]
      # polymers = frame[frame[:,1] == 1][:,2:]
      # polymers = polymers * (box[:,1] - box[:,0]) + box[:,0]
      polymers = polymers.reshape((nMon,-1,3),order='F')

      # Calculate COM and RG
      COM = polymers.mean(0)
      RG = (polymers - COM)**2
      RG = RG.sum(2)
      RG = np.sqrt(RG.mean(0))

      box = np.array([[-50,50],[-50,50]], dtype=float)
      hist = make2dHist(COM[:,:2], RG, box)

      # output
      if len(RG) != 0:
        # write COM, RG and surface atoms
        output = np.zeros((len(RG),4))
        output[:,:3] = COM
        output[:,3] = RG
        # writeChains(outFile, output, t, mode)
        writeChains(outFile, hist, t, mode)
        print('%d Chains written at timestep: %d' % (len(output),t))
      else:
        print 'No observed infiltration'
          
  print('Finished analyzing trajectory')

if __name__ == "__main__":  
  main()
