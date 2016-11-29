import numpy as np
from sys import argv, exit
import argparse
from conf_tools import readFrame, readTrj

NDIMS = 3
AREA = 100*100

def getArgs(argv):
  
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', '--trajFile')
  parser.add_argument('-o', '--output')
  parser.add_argument('-n', '--nFrames', type=int)
  parser.add_argument('-b', '--binSize', type=int)
  parser.add_argument('-r', type=float)

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
  above = data[data[:,2] > cutOff]
  below = data[data[:,2] <= cutOff]
 
  if below[-1,1] == cutOff:
    height = below[-1,0]
  else:
    height = interpolate(cutOff, below[-1,2::-2], above[0,2::-2])
  
  return height
 
def makeHistogram(atoms,R,binSize,limits):
  # make a histogram with different sized bins
  nBins = int((limits[1]-limits[0])/binSize)
  binSize = (limits[1]-limits[0])/nBins
  hist = np.zeros((nBins,3))
  hist[:,0] = np.arange(nBins)*binSize + limits[0]

  # calculate bin numbers for all atoms
  interior = np.logical_and(atoms[:,2] >= limits[0], atoms[:,2] < limits[1])
  interior = np.logical_and(interior, np.sqrt(atoms[:,0]**2 + atoms[:,1]**2) < R)
  atoms = atoms[interior]
  bins = ((atoms - limits[0]) / binSize).astype(int)
  # loop through bins and add atoms
  hist[:,1] = [(bins == bin).sum() for bin in range(nBins)]

  # normalize by volume
  hist[:,1] /= binSize*float(np.pi*R**2)
  
  # calculate cumulative density profile
  hist[:,2] = np.cumsum(hist[:,1]) / hist[:,1].sum()
  
  return hist
  
def main():
  
  # Command line args processing
  args = getArgs(argv)
  
  with open('density_'+args.output, 'w') as otp:
    otp.write('# Density profiles\n')
  with open('height_'+args.output, 'w') as otp:
    otp.write('# height profiles at two cut-offs\n')

  with open(args.trajFile,'r') as inp:
    for n in range(args.nFrames):
      ### .lammpstrj
      # read data and separate by atomtypes
      time, box, frame = readTrj(inp)
      frame = frame[frame[:,0].argsort()]
      #frame[:,2:] = frame[:,2:] * (box[:,1] - box[:,0]) + box[:,0]
      polymers = frame[frame[:,1] == 1][:,2:5]
      cylinder = frame[frame[:,1] == 3][:,2:]

      # ### .xyz
      # time, frame = readFrame(inp)
      # polymers = frame[frame[:,0] == 1][:,3]
      # cylinder = frame[frame[:,0] == 3][:,1:]

      #pmin = polymers.min()
      #polymers -= pmin
      #cylinder[:,2] -= pmin
      polymers[:,2] -= cylinder[:,2].min()
      cylinder[:,2] -= cylinder[:,2].min()

      # make box
      box = np.zeros((NDIMS,2))
      box[:,0] = cylinder.min(0)
      box[:,1] = cylinder.max(0)

      # calculate the density of the polymers in the packing
      density = makeHistogram(polymers, args.r, args.binSize, box[2])

      ## Integrate real density to the appropriate cut-off
      #cutOff = 0.85
      #height85 = integrateHeight(density, cutOff)

      ## Integrate real density to the appropriate cut-off
      #cutOff = 0.99
      #height99 = integrateHeight(density, cutOff)

      # output density and height values
      with open('density_'+args.output, 'a') as otp:
        header = '# time: %d\n#  z  density  cum_density' % time
        np.savetxt(otp, density, fmt='%.5f', header=header, footer='\n', comments='')
      #with open('height_'+args.output, 'a') as otp:
      #  np.savetxt(otp, [[time, height85, height99, density[1,0]]], fmt='%d %.5f %.5f %.5f',comments='')
      
  print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
