#!/usr/bin/python

import numpy as np
from sys import argv, exit
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from getopt import *
from matplotlib import pyplot
from pylab import *
from operator import lt, ge

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

def readFrame(file):
  line = file.readline().strip()
  nAtoms = int(line)
  atoms = np.zeros((nAtoms,4))

  line = file.readline()
  time = int(line.split()[-1])

  for i in range(nAtoms):
    line = file.readline().split()
    atoms[i,:] = np.array(line, dtype=float) 

  return time, atoms
 
def combineFrames(file, nCombine):
  
  time, atoms = readFrame(file)
  l = len(atoms)
  combined = np.zeros((l*nCombine,4))
  combined[0:l] = atoms
  for n in range(nCombine-1):
    _, atoms = readFrame(file)
    combined[l*(n+1):l*(n+2)] = atoms
    
  return time, combined
 
def centerAtoms(atoms):

  com = atoms.mean(0)
  atoms[:,:2] = atoms[:,:2] - com[:2]

  return atoms

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

  return interfacial # [interfacial[:,1] == rhoz[:,0][rhoz[:,0] < p[1]]]

def contactAngle(interface,box):
  # fit interfacial atom coordinates to a circle
  interface[:,1] -= box[2,0]
  R = lambda y0, z0: np.sqrt((interface[:,0] - y0)**2 + (interface[:,1] - z0)**2)
  # res = lambda z0: R(z0) - R(z0).mean()
  res = lambda args: R(*args) - R(*args).mean()
  popt, pcov = leastsq(res, [interface[:,0].mean(),-interface[:,1].mean()])
  # print popt
  # angle = np.arccos(popt[1]/popt[0])*180.0/np.pi
  angle = np.arccos(-popt[1]/R(*popt).mean())*180.0/np.pi

  return angle, popt[1], R(*popt).mean()

def plotInterface(atoms,box,interface,circle,t_angle):
  fig = pyplot.figure()
  ax = fig.add_subplot(1,1,1)
  # pyplot.plot(atoms[:,1],atoms[:,2],'o')
  # print interface
  pyplot.plot(interface[:,0],interface[:,1],'go')
  if circle is not None:
    fig.gca().add_artist(circle)
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

def plotRhoZ():
  pass
  # guess = [3,40,1]
  # p, cov = curve_fit(rho,rhoz[:,0],rhoz[:,1],guess)
  # pyplot.figure()
  # pyplot.plot(rhoz[:,0],rhoz[:,1],'o')
  # pyplot.plot(rhoz[:,0],rho(rhoz[:,0],*p))
  # pyplot.ylim(bottom=0)

def plotAngle(angles):
  angles = np.array(angles)
  angles = angles[angles[:,0] > 1000]
  pyplot.figure()
  pyplot.plot(angles[:,0],angles[:,1],'bo')
  pyplot.ylabel(r'Contact Angle, $\theta$ (deg)')
  pyplot.xlabel(r'time, t ($\tau_{LJ}$)')
  pyplot.title(r'$\epsilon$ = 0.50, $<\theta>$ = %.0f $\pm$ %.0f' % 
               (angles[:,1].mean(),angles[:,1].std()))

def writeData(filename, data):
  with open(filename, 'w') as file:
    file.write('# time angle\n')
    [file.write('%d %f\n' % tuple(line)) for line in data]

def main():
  
  trajFile, outFile, nSteps = getArgs(argv)

  nAve = 1
  circ = None
  frames = []
  with open(trajFile,'r') as inp:
    for n in range(nSteps):
      # How should block averaging of interfaces be done?
      # t, atoms = combineFrames(inp,nAve)
      t, atoms = readFrame(inp)

      # remove vaporized chains by uncommenting last two mask conditions
      mask = (atoms[:,0] == 1) # & (atoms[:,3] < 20) & (abs(atoms[:,2]) < 60)
      atoms = atoms[mask][:,1:]
      # Find the center of mass in the xy direction and center droplet
      atoms = centerAtoms(atoms)

      # make box
      box = np.array([atoms.min(0)-0.01,atoms.max(0)+0.01]).T
      box[1,:] += [-7,7]
      atoms[:,2] -= box[2,0]
      box[2,:] -= box[2,0]

      # find interfacial atoms
      interface = findCylinderInterface(atoms, box)
      # pyplot.show()
      # Calculate the contact angle
      if not ((n+1) % nAve):
        try:
          angle, z0, R = contactAngle(interface[(interface[:,1] < 30) & 
                                                (abs(interface[:,0]) < 100)], box)
        except RuntimeError:
          print 'No interface'
          circ = None
        except IndexError:
          print 'Empty interface'
        else:
          # pyplot.plot(interface[:,0],interface[:,1],'o')
          circ = pyplot.Circle((0,z0), radius=R, fill=False)
          frames.append([t*0.002,angle])

      #     plotInterface(atoms,box,interface,circ,(t,angle))
      # pyplot.show()

    plotAngle(frames)
    pyplot.show()
    writeData(outFile,frames)
 
    print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
