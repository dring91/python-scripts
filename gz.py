import numpy as np
from sys import argv, exit
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from getopt import *
from matplotlib import pyplot
from pylab import *
from operator import lt, ge
from conf_tools import readFrame
from pbc_tools import *

def getArgs(argv):
  
  try:
    opts, args = getopt(argv[1:],'t:o:n:')
  except GetoptError:
    print 'calcInterface.py -t <filename> -o <filename> -n <nSteps>'
    exit(2)
  for opt, arg in opts:
    if opt == '-t':
      trajFile = arg
    elif opt == '-o':
      confFile = arg
    elif opt == '-n':
      nSteps = arg

  return trajFile, confFile, int(nSteps)

def centerAtoms(atoms):

  com = atoms.mean(0)
  atoms[:,:2] = atoms[:,:2] - com[:2]

  return atoms, com[1]

def partitionAtoms(atoms, box, nSlices):

  dz = box[2,1]/nSlices
  z = dz * (np.arange(nSlices) + 0.5)
  slices = [[] for i in range(nSlices)]
  for atom in atoms:
    slice = int(atom[2]/dz)
    slices[slice].append(atom[1])

  return slices,z,dz

def createHistogram(slice, box, nBins, dz):

  dy = (box[1,1]-box[1,0])/nBins
  density = np.zeros((nBins, 2))
  density[:,0] = dy * (np.arange(nBins) + 0.5) + box[1,0]
  for atom in slice:
    index = int((atom-box[1,0])/dy)
    density[index,1] += 1.0
  density[:,1] = density[:,1]/dy/dz/(box[0,1]-box[0,0]) 

  return density
  
def rho(z, dl, ze, l):
    
  return dl/2*(1-np.tanh(2*(z-ze)/l))
 
def fitDensity(op, z, density):
  
  if op(1,2):
    guess = [0.8,-60,-1]
    fstr = 'bo'
  else:
    guess = [0.8,60,1]
    fstr = 'go'
  mask = op(density[:,0],0)
  if sum(mask) < 3:
    return None

  try:
    params, cov = curve_fit(rho,density[mask][:,0],density[mask][:,1],guess)
  except RuntimeError:
    print 'z = %.2f slice not fit' % z
  else:
    # Save the "half-density" value to its z-slice
    if 0.7 < params[0] < 1: # cov[1,1] < 10**3 and not np.isinf(cov[1,1]):      
      # plotRhoY(z, density,mask,fstr,params)
      return params[1]
    
def findCylinderInterface(atoms, box):

  # define resolution of the z (slices) and y (bins) coordinates
  # There are ~60k atoms per frame
  # 60,000 / 15 / 50 = 80 atoms per cell on average
  # for 5 combined frames, the value is 400
  zRes = 15
  yRes = 50
  # sort atoms according to slice (ragged list)
  slices, zs, dz = partitionAtoms(atoms, box, zRes)

  interfacial = []
  rhoz = []
  # iterate over z slices
  for index, (slice, z) in enumerate(zip(slices,zs)):
    # for the nth slice, calculate the radial density
    density = createHistogram(slice, box, yRes, dz)
    # fit density to sigmoidal function
    yl = fitDensity(lt,z,density)
    if yl is not None:
      interfacial.append([yl,z])
    yr = fitDensity(ge,z,density)
    if yr is not None:
      interfacial.append([yr,z])
    # if yl is not None and yr is not None:
      # rhoz.append([z,density[(yl < density[:,0]) & (density[:,0] < yr)][:,1].mean()])
    # rhoz.append([z,density[density[:,1] > 0][:,1].mean()])
  # rhoz = np.array(rhoz)
  interfacial = np.array(interfacial)

  return interfacial

def plotInterface(atoms,box,interface):
  fig = pyplot.figure()
  ax = fig.add_subplot(1,1,1)
  # pyplot.plot(atoms[:,1],atoms[:,2],'o')
  # print interface
  pyplot.plot(interface[:,0],interface[:,1],'go')
  pyplot.ylabel(r'Height, z ($\sigma$)')
  pyplot.xlabel(r'Width, y ($\sigma$)')
  pyplot.ylim(box[2,:])
  pyplot.xlim(box[1,:])
  pyplot.title(r't = %d, $\theta$ = %f' % t_angle)
  ax.set_aspect('equal','datalim')
 
def plotRhoY(z,density,mask,fstr,params):
  # pyplot.figure()
  pyplot.title(r'z = %f ($\sigma$)' % z)
  pyplot.plot(density[mask][:,0],density[mask][:,1],fstr)
  pyplot.xlim((-60,60))
  pyplot.xlabel(r'y coordinate ($\sigma$)')
  pyplot.ylabel(r'Density, $\rho$(y) ($n/\sigma^{3}$)')
  pyplot.plot(density[mask][:,0],rho(density[mask][:,0],*params)) 

def plotRhoZ(rhoz):
  # guess = [3,40,1]
  # p, cov = curve_fit(rho,rhoz[:,0],rhoz[:,1],guess)
  pyplot.figure()
  pyplot.plot(rhoz[:,0],rhoz[:,1],'o')
  # pyplot.plot(rhoz[:,0],rho(rhoz[:,0],*p))
  pyplot.ylim(bottom=0)

def writeData(filename, data):
  with open(filename, 'w') as file:
    file.write('# time angle\n')
    [file.write('%d %f\n' % tuple(line)) for line in data]

def calcPCF(surf, poly, box, nBins):
  rc = box[2,1]
  rho = 0.8
  binSize = rc / nBins
  pcf = np.zeros((nBins,2))
  pcf[:,0] = (np.arange(nBins)+0.5)*binSize

  for s in surf:
    dr = poly-s
    dr[:,:2] = PBC(dr[:,:2],dr[:,:2],box[:2,:])
    r = np.sqrt(dr[:,0]**2+dr[:,1]**2+dr[:,2]**2)
    mask = (r < rc)
    bins = (r[mask] / binSize).astype(int)
    pcf[bins,1] += 1
  pcf[1:,1] /= len(surf)*rho*4.0*np.pi*(pcf[1:,0]**2 - pcf[:nBins-1,0]**2)*binSize

  return pcf

def main():
  
  trajFile, outFile, nSteps = getArgs(argv)

  nBins = 180
  pcf = np.zeros((nBins,2))
  with open(trajFile,'r') as inp:
    for n in range(nSteps):
      for k in range(50):
        t, frame = readFrame(inp)

      mask = lambda i: (frame[:,0] == i)
      poly = frame[mask(1)][:,1:]
      surf = frame[mask(2)][:,1:]
      # Find the center of mass in the xy direction and center droplet
      # poly, center = centerAtoms(poly)

      # make box
      box = np.ones((3,2))
      box[:,0] = frame[:,1:].min(0) #-0.01
      box[:,1] = frame[:,1:].max(0) #+0.01
      # box[1,:] += [-7,7]

      # atoms[:,2] -= box[2,0]
      # box[2,:] -= box[2,0]

      # find interfacial atoms
      # interface = findCylinderInterface(atoms, box)

      # calculate volume per horizontal 'slice'
      # calculate pair correlation function
      pcf += calcPCF(surf, poly, box, nBins)
    plt.plot(pcf[:,0], pcf[:,1])
    plt.show()

    print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
