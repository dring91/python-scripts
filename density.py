import numpy as np
from sys import argv, exit
from argparse import ArgumentParser
from conf_tools import readFrame, readTrj
from scipy.optimize import curve_fit, OptimizeWarning
import matplotlib.pyplot as plt
from functools import partial
from Jacobian import jac

NDIMS = 3
AREA = 100*100

def makeHistogram(atoms,R,binSize,limits):
  # make a histogram with different sized bins
  nBins = int((limits[1]-limits[0])/binSize)
  binSize = (limits[1]-limits[0])/nBins
  hist = np.zeros((nBins,4))
  hist[:,0] = (np.arange(nBins)+0.5)*binSize + limits[0]

  # calculate bin numbers for all atoms
  interior = np.logical_and(atoms[:,2] >= limits[0], atoms[:,2] < limits[1])
  interior = np.logical_and(interior, np.sqrt(atoms[:,0]**2 + atoms[:,1]**2) < R)
  atoms = atoms[interior]
  bins = ((atoms[:,2] - limits[0]) / binSize).astype(int)
  # loop through bins and add atoms
  hist[:,1] = [(bins == bin).sum() for bin in range(nBins)]

  # normalize by volume
  hist[:,2] = hist[:,1] / (binSize*float(np.pi*(R-0.5)**2))
  
  # calculate cumulative density profile
  hist[:,3] = np.cumsum(hist[:,2]) / hist[:,2].sum()
  
  return hist
  
def sigmoidal(z, rhol, rhov, z0, d):
  return 0.5 * (rhol + rhov) - 0.5 * (rhol - rhov) * np.tanh(2*(z - z0)/d)

no_vapor = partial(sigmoidal, rhov=0) 
const_interface_width = partial(no_vapor, d=2) 
const_liquid_density = partial(const_interface_width, rhol=0.8) 

def main():

  # Command line args processing
  parser = ArgumentParser()
  parser.add_argument('-i', '--input')
  parser.add_argument('-o', '--output')
  parser.add_argument('-n', '--nFrames', type=int)
  parser.add_argument('-b', '--binSize', type=float, 
                      help='estimate of the histogram bin size')
  parser.add_argument('-r', type=float, help='radius of cylinder for calculating density')

  args = parser.parse_args()
  
  outFile = '%s_B%d' % (args.output, args.binSize)
  
  with open('density_'+outFile, 'w') as otp:
    otp.write('# Density profiles\n')
  with open('height_'+outFile, 'w') as otp:
    otp.write('# height profiles at two cut-offs\n')
    otp.write('time rhol rhov height interface polymers height85 height99\n')
  with open('errors_'+outFile, 'w') as otp:
    otp.write('# errors from fit\n')

  height = 0
  with open(args.input,'r') as inp:
    for n in range(args.nFrames):
      ### .lammpstrj
      # read data and separate by atomtypes
      time, box, frame = readTrj(inp)
      frame = frame[frame[:,0].argsort()]
      #frame[:,2:] = frame[:,2:] * (box[:,1] - box[:,0]) + box[:,0]
      polymers = frame[frame[:,1] == 2][:,2:5]
      cylinder = frame[frame[:,1] == 4][:,2:]

      # ### .xyz
      # time, frame = readFrame(inp)
      # polymers = frame[frame[:,0] == 1][:,3]
      # cylinder = frame[frame[:,0] == 3][:,1:]

      # vertical zero at the bottom of the cylinder
      #polymers[:,2] -= cylinder[:,2].min()
      #cylinder[:,2] -= cylinder[:,2].min()

      # make box
      box = np.zeros((NDIMS,2))
      box[:,0] = cylinder.min(0)
      box[:,1] = cylinder.max(0)

      # calculate the density of the polymers in the packing
      density = makeHistogram(polymers, args.r, args.binSize, box[2])
      density[:,0] -= cylinder[:,2].min()

      # Calculate fluid height via sigmoidal fitting
      guess = [0.8, 0.05, height, 5]
      bounds = (0, 105)

      # trf method
      params, covariance = curve_fit(sigmoidal, density[:,0], density[:,2], guess, bounds=bounds)
      guess = [0.8, height, 5]
      params, covariance = curve_fit(no_vapor, density[:,0], density[:,2], guess, bounds=bounds)
      guess = [0.8, height]
      params, covariance = curve_fit(const_interface_width, density[:,0], density[:,2], guess, bounds=bounds)
      guess = [height]
      params, covariance = curve_fit(const_liquid_density, density[:,0], density[:,2], guess, bounds=bounds)

      # dogbox method
      params, covariance = curve_fit(sigmoidal, density[:,0], density[:,2], guess, bounds=bounds, method='dogbox')
      guess = [0.8, height, 5]
      params, covariance = curve_fit(no_vapor, density[:,0], density[:,2], guess, bounds=bounds, method='dogbox')
      guess = [0.8, height]
      params, covariance = curve_fit(const_interface_width, density[:,0], density[:,2], guess, bounds=bounds, method='dogbox')
      guess = [height]
      params, covariance = curve_fit(const_liquid_density, density[:,0], density[:,2], guess, bounds=bounds, method='dogbox')

      # lm method
      params, covariance = curve_fit(sigmoidal, density[:,0], density[:,2], guess, method='lm')
      guess = [0.8, height, 5]
      params, covariance = curve_fit(no_vapor, density[:,0], density[:,2], guess, method='lm')
      guess = [0.8, height]
      params, covariance = curve_fit(const_interface_width, density[:,0], density[:,2], guess, method='lm')
      guess = [height]
      params, covariance = curve_fit(const_liquid_density, density[:,0], density[:,2], guess, method='lm')

      height = params[2]

      # MSE (Mean Square Error) is a measure of fit
      #MSE = np.einsum('i...,...i',covariance,jac())
      #MSE = np.einsum('i...,i...',jac(),MSE)

      # Interpolate integrated density to the appropriate cut-off
      cutOffs = [0.85, 0.99] 
      heights = np.interp(cutOffs, density[:,3], density[:,0])

      # output density and height values
      with open('density_'+outFile, 'a') as otp:
        header = '# time: %d\n#  z  counts  density  cum_density' % time
        np.savetxt(otp, density, fmt='%.5f', header=header, footer='\n', comments='')
      with open('height_'+outFile, 'a') as otp:
        np.savetxt(otp, [[time] + params.tolist() + [density[1,0]] + heights.tolist()], 
                   fmt='%.5f')
      with open('errors_'+outFile, 'a') as otp:
        try:
          np.savetxt(otp, covariance.reshape((1,len(params)**2)))
        except AttributeError:
          np.savetxt(otp, np.zeros((1,len(params)**2)))
      
  print('Finished analyzing trajectory')

if __name__ == '__main__':
  main()
