import numpy as np
from sys import argv, exit
import argparse
from conf_tools import readFrame, readTrj
from scipy.optimize import curve_fit, OptimizeWarning

NDIMS = 3
AREA = 100*100

def getArgs(argv):
  
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', '--trajFile')
  parser.add_argument('-o', '--output')
  parser.add_argument('-n', '--nFrames', type=int)
  parser.add_argument('-b', '--binSize', type=float)
  parser.add_argument('-r', type=float)

  return parser.parse_args()
  
def interpolate(x, low, high):
  if x <= high[0] and x >= low[0]:
    slope = (high[1] - low[1])/(high[0] - low[0])
    return slope*(x - low[0]) + low[1]
  else:
    print 'x out of range'
    exit(0)

def integrateHeight(data, cutOff, col=2):
  # filter values above and below cutoff
  above = data[data[:,col] > cutOff]
  below = data[data[:,col] <= cutOff]
 
  try:
    height = interpolate(cutOff, below[-1,col::-col], above[0,col::-col])
  except IndexError:
    height = 0
  
  return height
 
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

  # params = [rhol, rhov, z0, d]

  return 0.5 * (rhol + rhov) - 0.5 * (rhol - rhov) * np.tanh(2*(z - z0)/d)

def jac(z, rhol, rhov, z0, d):

  sech = lambda x: 1/np.cosh(x)

  J = [0.5 - 0.5 * np.tanh(2*(z-z0)/d),
       0.5 * np.tanh(2*(z-z0)/d) - 0.5,
       (rhol - rhov)*z0/d*sech(2*(z-z0)/d)**2,
       (rhol - rhov)/d**2*(z-z0)*sech(2*(z-z0)/d)**2]

  return J

def main():
  
  # Command line args processing
  args = getArgs(argv)
  outFile = '%s_b%.1f' % (args.output, args.binSize)
  
  with open('density_'+outFile, 'w') as otp:
    otp.write('# Density profiles\n')
  with open('height_'+outFile, 'w') as otp:
    otp.write('# height profiles at two cut-offs\n')
  with open('cut-offs_'+outFile, 'w') as otp:
    otp.write('# heights calculated for multiple cut-offs\n')

  height = 0
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
      guess = [0.8, 0.1, height, 2]
      bounds = [0, 105]

      try:
        params, errors = curve_fit(sigmoidal, density[:,0], density[:,2], guess, bounds=bounds, method='trf')
      except RuntimeError:
        params, errors = curve_fit(sigmoidal, density[:,0], density[:,2], guess, method='lm')

      height = params[2]

      # Integrate real density to the appropriate cut-off
      cutOff = 0.85
      height85 = integrateHeight(density, cutOff, 3)

      # Integrate real density to the appropriate cut-off
      cutOff = 0.99
      height99 = integrateHeight(density, cutOff, 3)

      # compute different cut-off values
      cutOffs = np.arange(0.85,0.99,0.01)
      heights = [integrateHeight(density, cutOff, 3) for cutOff in cutOffs]

      # output density and height values
      with open('density_'+outFile, 'a') as otp:
        header = '# time: %d\n#  z  counts  density  cum_density' % time
        np.savetxt(otp, density, fmt='%.5f', header=header, footer='\n', comments='')
      with open('height_'+outFile, 'a') as otp:
        np.savetxt(otp, [[time] + params.tolist() + [height85, height99, density[1,0]]], 
                        fmt='%d %.5f %.5f %.5f %.5f %.5f %.5f %.5f',
                        header='time rhol rhov height interface height85 height99 polymers')
      with open('cut-offs_'+outFile, 'a') as otp:
        np.savetxt(otp, [[time] + heights + params.tolist()], fmt='%.5f',
                        header='time rhol rhov height interface height85 height99 polymers')
      
  print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
