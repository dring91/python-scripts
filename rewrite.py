import numpy as np
from sys import argv, exit
from getopt import *
import matplotlib.pyplot as plt
from pbc_tools import *
from conf_tools import readFrame, readTrj

def getArgs(argv):
  
  try:
    opts, args = getopt(argv[1:],'t:o:n:')
  except GetoptError:
    print 'randomWalk.py -t <filename>'
    exit(2)
  for opt, arg in opts:
    if opt == '-t':
      trajFile = arg
    if opt == '-o':
      outFile = arg
    if opt == '-n':
      nFrames = arg

  return trajFile, outFile, int(nFrames)
 
def main():

  trajFile, outFile, nFrames = getArgs(argv)
  atomtype = 1
  particles = 3
  nSkip = 1 # nFrames / 100
  nMon = 5

  with open(outFile+'_unwrapped.xyz', 'w') as otp:
    otp.write('')

  with open(trajFile,'r') as inp:
    for n in range(nFrames/nSkip):
      ### lammmpstrj ###
      # # read and filter data
      # for s in range(nSkip):
      #   time, box, frame = readTrj(inp)
      # zero = frame[frame[:,1] == particles][:,4].min()
      # frame = frame[np.argsort(frame[:,0])]
      # frame = frame[frame[:,1] == atomtype][:,2:]
      # Lxy = box[:,1] - box[:,0]
      # frame = frame * Lxy + box[:,0]

      ### xyz ###
      # read and filter data
      for s in range(nSkip):
        time, frame = readFrame(inp)
      zero = frame[frame[:,0] == particles][:,3].min()
      frame = frame[frame[:,0] == atomtype][:,1:]
      frame[:,2] -= zero

      # build box
      box = np.zeros((2,2))
      box[:,0] = frame[:,:2].min(0)
      box[:,1] = frame[:,:2].max(0)
      Lxy = box[:,1] - box[:,0]
      frame = unwrap(frame.reshape((-1,nMon,3)), box).reshape((-1,3))
      if n == 0:
        periods = np.zeros_like(frame[:,:2])
        move = np.zeros_like(frame[:,:2])
      elif n > 0:
        move = frame[:,:2] + periods * Lxy[:2] - oldFrame[:,:2]
        periods += (move < box[:2,0])
        periods -= (move >= box[:2,1])
        frame[:,:2] += periods * Lxy[:2]
      # if differences are greater than L/2 in any frame, add 1
      oldFrame = np.copy(frame)

      with open(outFile+'_unwrapped.xyz', 'a') as otp:
        otp.write('%d\n' % len(frame))
        otp.write('Atoms. Timestep: %d\n' % time)
        np.savetxt(otp,frame,'1 %.5f %.5f %.5f')
     
    print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
