import numpy as np
from scipy import ndimage
from scipy.optimize import curve_fit
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
  parser.add_argument('-b', '--binSize', nargs=2, type=float)
  parser.add_argument('-r', type=float)

  return parser.parse_args()
  
def makeHistogram(atoms,R,binSize,limits):
  # make a histogram with different sized bins
  L = limits[1]-limits[0]
  nBins = (int(L/binSize[0]), int(R/binSize[1]))
  binSize = (L/nBins[0], R/nBins[1])
  hist = np.zeros((nBins[0]+1, nBins[1]+1))
  hist[1:,0] = (np.arange(nBins[0])+0.5)*binSize[0]+limits[0]
  hist[0,1:] = (np.arange(nBins[1])+0.5)*binSize[1]

  # calculate bin numbers for all atoms
  interior = np.logical_and(atoms[:,2] >= limits[0], atoms[:,2] < limits[1])
  interior = np.logical_and(interior, np.sqrt(atoms[:,0]**2 + atoms[:,1]**2) < R)
  atoms = atoms[interior]
  bins = np.zeros((len(atoms),2))
  bins[:,0] = (atoms[:,2] - limits[0]) / binSize[0]
  bins[:,1] = np.sqrt(atoms[:,0]**2 + atoms[:,1]**2) / binSize[1]
  bins = bins.astype(int)

  # loop through bins and add atoms
  hist[1:,1:] = [[np.logical_and(bins[:,0] == l,bins[:,1] == r).sum() for r in range(nBins[1])] for l in range(nBins[0])]

  # normalize by volume
  hist[1:,1:] = hist[1:,1:] / (np.pi*binSize[0]*binSize[1]**2*(2*np.arange(nBins[1])+1))
  
  return hist

def sigmoidal(z, rhol, rhov, z0, d):

  return 0.5 * (rhol + rhov) - 0.5 * (rhol - rhov) * np.tanh(2*(z - z0)/d)

def jac(z, rhol, rhov, z0, d):

  return 0.5 * (rhol - rhov0) * np.sech(2*(z - z0)/d)**2 * (2*z/d)
  
def main():
  
  # Command line args processing
  args = getArgs(argv)
  
  with open('density_'+args.output, 'w') as otp:
    otp.write('# Density profiles\n')

  with open(args.trajFile,'r') as inp:
    for n in range(args.nFrames):
      ### .lammpstrj
      # read data and separate by atomtypes
      time, box, frame = readTrj(inp)
      frame = frame[frame[:,0].argsort()]
      polymers = frame[frame[:,1] == 1][:,2:5]
      cylinder = frame[frame[:,1] == 3][:,2:5]

      # remove center of mass
      com = cylinder[:,:2].mean(0)
      cylinder[:,:2] -= com

      # make box
      box = np.zeros((NDIMS,2))
      box[:,0] = cylinder.min(0)
      box[:,1] = cylinder.max(0)

      # calculate the density of the polymers in the packing
      density = makeHistogram(polymers, args.r, args.binSize, box[2])

      # calculate interface from sigmoidal function
      guess = [0.8, 0.1, 50, 2]
      interface = []
      for data in density.T:
        params, errors = curve_fit(sigmoidal, density[1:,0], data[1:], guess)
        interface.append(params)
      interface = np.array(interface) 

      # calculate meniscus from Sobel filter
      #sx = ndimage.sobel(density[1:,1:], axis=1, mode='reflect')
      #sy = ndimage.sobel(density[1:,1:], axis=0, mode='reflect')
      #meniscus = np.hypot(sx, sy)
      #edge = np.argwhere(meniscus > 1.5)

      # output density and height values
      with open('density_'+args.output, 'a') as otp:
        header = '# time: %d\n#  z  counts  density' % time
        #np.savetxt(otp, density[1:,1:], fmt='%.5f', header=header, footer='\n', comments='')
        #np.savetxt(otp, meniscus, fmt='%.5f', header=header, footer='\n', comments='')
        np.savetxt(otp, interface, fmt='%.5f', header=header, footer='\n', comments='')
     
  print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
