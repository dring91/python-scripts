import numpy as np
from conf_tools import readTrj,iTrj
import argparse

def main():
  """ A script to calculate the rouse modes and their contribution to viscosity """

  parser = argparse.ArgumentParser()
  parser.add_argument('--input','-i')
  parser.add_argument('--output','-o')
  parser.add_argument('--frames','-n',type=int)
  parser.add_argument('--length','-l',type=int)
  args = parser.parse_args()

  with open(args.input, 'r') as file:
    #for time, box, frame in iTrj(file):
    for n in xrange(args.frames):
      time, box, frame = readTrj(file)

      frame = frame[np.argsort(frame[:,0])]
      polymers = frame[frame[:,1] == 1][:,2:]
      nChains = polymers.shape[0]/args.length
      polymers = polymers.reshape((nChains,args.length,3))

      stress = 0
      for p in xrange(1,args.length):
        modes = np.cos((np.arange(args.length)+0.5)*p*np.pi/args.length)
        modes = np.einsum('ijk,j->ik',polymers,modes) / args.length
        stress += np.einsum('i,i',modes[:,0],modes[:,1])
      print stress

if __name__ == '__main__':
  main()
