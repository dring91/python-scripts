import numpy as np
from sys import argv, exit
import argparse
from conf_tools import readFrame, readTrj
from insert_tools import *

NDIMS = 3
AREA = 100*100
RC = 0

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
  bins = ((polymers - box[2,0]) / dSlice).astype(int)
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

def getArgs(argv):
  
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', '--trajFile')
  parser.add_argument('-o', '--output')
  parser.add_argument('-n', '--nFrames', type=int)
  parser.add_argument('-i', '--nInsert', type=int)
  parser.add_argument('-s', '--nSlices', type=int)

  return parser.parse_args()
  
def interpolate(x, low, high):
  if x <= high[0] and x >= low[0]:
    slope = (high[1] - low[1])/(high[0] - low[0])
    return slope*(x - low[0]) + low[1]
  else:
    print 'x out of range'
    exit(0)

def integrateHeight(data, cutOff):
  # filter values above and below cutoff
  above = data[data[:,3] > cutOff]
  below = data[data[:,3] <= cutOff]
 
  if below[-1,1] == cutOff:
    height = below[-1,0]
  else:
    height = interpolate(cutOff, below[-1,3::-3], above[0,3::-3])
  
  return height
 
def makeHistogram(porosity,atoms,nBins,limits):
  # make a histogram with different sized bins
  hist = np.zeros((nBins+2,4))
  binSize = (limits[1]-limits[0])/float(nBins)
  hist[1:,0] = np.arange(nBins+1)*binSize + limits[0]
  hist[1,1] = (atoms < limits[0]).sum()

  # calculate bin numbers for all atoms
  atoms = atoms[np.logical_and(atoms >= limits[0], atoms < limits[1])]
  bins = ((atoms - limits[0]) / binSize).astype(int)
  # loop through bins and add atoms
  hist[2:,1] = [(bins == bin).sum() for bin in range(2,nBins+2)]

  # normalize by volume
  hist[:,1] /= float(AREA)
  hist[2:,1] /= binSize
  hist[1,1] /= limits[0]
  hist[0,1] = hist[1,1]
  
  # divide by the porosity to recovery the real density
  hist[2:,1] /= porosity[:,1]
  hist[2:,2] = porosity[:,1]

  # calculate cumulative density profile
  widths = hist[1:,0]-hist[:-1,0]
  hist[1:,3] = np.cumsum(hist[1:,1]*widths) / (hist[1:,1]*widths).sum()
  
  return hist
  
def main():
  
  # Command line args processing
  args = getArgs(argv)
  
  outFile = '%s_I%d_B%d' % (args.output, args.nInsert, args.nSlices)

  # params = (a,b,r_cut,nAtomsPerPart,nPart)
  #params = (25*0.5, 50*0.5, RC, 4684, 54)
  params = {'a':25*0.5, 'b':50*0.5, 'r':1}
  ns = {'nPart':54, 'nAtoms':4684, 'nDim':3, 'nBins':args.nSlices, 'nInsert':args.nInsert}

  with open('density_'+outFile, 'w') as otp:
    otp.write('# Density profiles\n')
  with open('height_'+outFile, 'w') as otp:
    otp.write('# height profiles at two cut-offs\n')

  with open(args.trajFile,'r') as inp:
    for n in range(args.nFrames):
      ### .lammpstrj
      # read data and separate by atomtypes
      time, box, frame = readTrj(inp)
      frame = frame[frame[:,0].argsort()]
      frame[:,2:] = frame[:,2:] * (box[:,1] - box[:,0]) + box[:,0]
      polymers = frame[frame[:,1] == 1][:,4]
      particles = frame[frame[:,1] >= 3][:,2:]

      # ### .xyz
      # time, frame = readFrame(inp)
      # polymers = frame[frame[:,0] == 1][:,3]
      # particles = frame[frame[:,0] >= 3][:,1:]

      #pmin = polymers.min()
      #polymers -= pmin
      #particles[:,2] -= pmin

      # make box
      box = np.zeros((NDIMS,2))
      box[:,0] = particles.min(0)-0.5
      box[:,1] = particles.max(0)+0.5

      # Calculate the packing porosity
      porosity = calc_porosity(particles, polymers, box, ns, params)
      #porosity = calcPorosity(particles, box, args.nInsert, args.nSlices, params)
      # calculate the superficial density of the polymers in the packing
      density = makeHistogram(porosity, polymers, args.nSlices, box[2])

      # # Integrate real density to the appropriate cut-off
      # cutOff = 0.85
      # height85 = integrateHeight(density, cutOff)

      # # Integrate real density to the appropriate cut-off
      # cutOff = 0.99
      # height99 = integrateHeight(density, cutOff)

      # output density and height values
      with open('density_'+outFile, 'a') as otp:
        header = '# time: %d\n#  z  density  porosity' % time
        np.savetxt(otp, density, fmt='%.5f', header=header, footer='\n', comments='')
#      with open('height_'+outFile, 'a') as otp:
#        np.savetxt(otp, [[time, height85, height99, density[1,0]]], fmt='%d %.5f %.5f %.5f',comments='')
      
  print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
