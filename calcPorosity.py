import numpy as np
from sys import argv, exit
from future_builtins import zip
import argparse

from conf_tools import readTrj
from pbc_tools import *
from insert_tools import *

def getArgs(argv):

  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input")
  parser.add_argument("-o", "--output")
  parser.add_argument("-f", "--frames", type=int)
  parser.add_argument("-p", "--particles", type=int)
  parser.add_argument("-b", "--bins", type=int)
  
  return parser.parse_args()

def calc_porosity(particles, polymers, box, ns, params):
  # calculate bin size
  particles = particles.reshape((ns['nPart'],ns['nAtoms'],ns['nDim']))
  dSlice = (box[2,1] - box[2,0]) / ns['nBins']

  # initialize porosity
  porosity = np.zeros((ns['nBins'],3))
  porosity[:,0] = (np.arange(ns['nBins']) + 0.5) * dSlice + box[2,0]

  # unwrap particle coordinates
  particles = unwrap(particles, box)

  # find the transforms to translate and rotate particles
  rotations, centers = calcTransformations(particles, box)
  
  # Select insertion points
  points = np.random.rand(ns['nInsert'],3)
  points = points*(box[:,1] - box[:,0]) + box[:,0]
  
  # test for overlap with nanoparticles
  test = np.zeros_like(points[:,0],dtype=bool)
  for c,r in zip(centers, rotations):
    test = insert(points, c, r, box, test, params)

  # calculate the superficial polymer density
  volume = dSlice * (box[0,1] - box[0,0]) * (box[1,1] - box[1,0])
  bead = np.pi/6
  bins = ((polymers[:,2] - box[2,0]) / dSlice).astype(int)
  superficial = np.array([(bins == bin).sum()*bead/volume for bin in range(ns['nBins'])])

  # partition points by vertical height
  bins = ((points[:,2] - box[2,0]) / dSlice).astype(int)
  #porosity[:,1] = np.zeros(ns['nBins'])

  # calculate the ratio of accepted to total insertions
  for bin in range(ns['nBins']):
    insertions = np.logical_not(test[bins == bin])
    porosity[bin,1] = insertions.sum() / float(len(insertions))

  porosity[:,2] = porosity[:,1] - superficial

  return porosity

def main():
  
  # Command line args processing
  args = getArgs(argv)
  
  outFile = '%s_%d_%d' % (args.output, args.particles, args.bins)

  params = {'a':25*0.5, 'b':50*0.5, 'r':1}
  ns = {'nPart':54, 'nAtoms':4684, 'nDim':3, 'nBins':args.bins, 'nInsert':args.particles}

  with open(outFile, 'w') as otp:
    otp.write('# Porosity for N = 5 oversaturated infiltration\n')

  with open(args.input,'r') as inp:
    for n in range(args.frames):
      # read and filter data from xyz file
      #time, frame = readFrame(inp)
      #frame = frame[frame[:,0] >= 3][:,1:]

      # read and filter data from traj file
      time, box, frame = readTrj(inp)
      frame[:,2:] = frame[:,2:]*(box[:,1] - box[:,0]) + box[:,0]
      frame = frame[np.argsort(frame[:,0])]
      particles = frame[frame[:,1] == 3][:,2:]
      polymers = frame[frame[:,1] == 1][:,2:]

      # make box
      box = np.zeros((3,2))
      box[:,0] = particles.min(0)-0.5
      box[:,1] = particles.max(0)+0.5

      porosity = calc_porosity(particles, polymers, box, ns, params)
      
      with open(outFile, 'a') as otp:
        np.savetxt(otp, 
                   porosity,
                   fmt='%.5f',
                   header='# time: %d\n#   z (sigma)   porosity\n' % time,
                   footer='\n',
                   comments='')
      
    print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
