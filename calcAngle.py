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
    print 'randomWalk.py -t <filename> -o <filename> -n <chainLength>'
    exit(2)
  for opt, arg in opts:
    if opt == '-t':
      trajFile = arg
    elif opt == '-o':
      confFile = arg
    elif opt == '-n':
      chainLength = arg

  return trajFile, confFile, int(chainLength)

def readTraj(file, atype, time):
  
  traj = {}
  for line in file:
    L = line.split()
    if len(L) > 1 and L[1] == 'Timestep:':
      t = int(L[2])
      if t > time:
        break
      traj[t] = []
    if len(L) > 0 and L[0] in set(atype):
      traj[t].append(L[1:])
    
  return traj
  
def combineFrames(traj,nAve):
  
  times = traj.keys()
  combined = {}
  tc = 0
  for t in sorted(times):
    if not t % nAve:
      tc = t
      combined[tc] = []
    combined[tc].extend(traj[t])
    
  return combined

def centerAtoms(atoms):

  com = atoms.mean(0)
  atoms[:,:2] = atoms[:,:2] - com[:2]

  return atoms

def partitionAtoms(atoms, box, nSlices):

  dy = (box[2,1]-box[2,0])/nSlices
  y = dy * (np.arange(nSlices) + 0.5) + box[2,0]
  slices = [[] for i in range(nSlices)]
  for atom in atoms:
    slice = int((atom[2]-box[2,0])/dy)
    slices[slice].append(atom[1])

  return slices,y,dy

def createHistogram(slice, box, nBins, dz):

  dr = (box[1,1]-box[1,0])/nBins
  density = np.zeros((nBins, 2))
  density[:,0] = dr * (np.arange(nBins) + 0.5) + box[1,0]
  for atom in slice:
    index = int((atom-box[1,0])/dr)
    density[index,1] += 1.0
  density[:,1] = density[:,1]/dr/dz/(box[0,1]-box[0,0]) 

  return density
  
def rho(z, dl, ze, l):
    
  return dl/2*(1-np.tanh(2*(z-ze)/l))
  
def grho(z, dl, ze, l):
  gradient = [rho(z,dl,ze,l)/dl, dl/l*(1-np.tanh(2*(z-ze)/l)**2), dl*(z-ze)/(l**2)*(1-np.tanh(2*(z-ze)/l)**2)]
  return gradient
  
def fitDensity(op, density):
  
  if op(1,2):
    guess = [0.7,-40,-1]
    fstr = 'bo'
  else:
    guess = [0.7,40,1]
    fstr = 'go'
  mask = op(density[:,0],0)
  if sum(mask) < 3:
    return None

  # pyplot.plot(density[mask][:,0],density[mask][:,1],fstr)
  # pyplot.xlim((-60,60))
  # pyplot.xlabel(r'y coordinate ($\sigma$)')
  # pyplot.ylabel(r'Density, $\rho$(y) ($n/\sigma^{3}$)')
    
  try:
    popt, pcov = curve_fit(rho,density[mask][:,0],density[mask][:,1],guess)
  except RuntimeError:
    print 'slice not fit'
  else:
    # Save the "half-density" value to its z-slice
    if pcov[1,1] < 10**3 and not np.isinf(pcov[1,1]):      
      # grad = np.transpose(grho(density[mask][:,0],popt[0],popt[1],popt[2]))
      # fiterr = np.inner(np.transpose(np.inner(pcov,grad)),grad)
      # pyplot.plot(density[mask][:,0],rho(density[mask][:,0],popt[0],popt[1],popt[2]))
      # pyplot.errorbar(density[mask][:,0],rho(density[mask][:,0],popt[0],popt[1],popt[2]),yerr=np.sqrt(pcov[1,1]))
      return popt[1]
    
def findCylinderInterface(atoms, box):
  # Find the center of mass in the xy direction and center droplet
  atoms = centerAtoms(atoms)

  # define resolution of the z (slices) and y (bins) coordinates
  zRes = 15
  yRes = 50
  # sort atoms according to slice (ragged list)
  slices, zs, dz = partitionAtoms(atoms, box, zRes)

  interfacial = []
  # iterate over z slices
  for index, (slice, z) in enumerate(zip(slices,zs)):
    # for the nth slice, calculate the radial density
    density = createHistogram(slice, box, yRes, dz)
    # pyplot.plot(density[:,0], density[:,1], 'o')
    # fit density to sigmoidal function
    # pyplot.figure(index)
    # pyplot.title(r'z = %f ($\sigma$)' % (z-box[2,0]))
    y = fitDensity(lt,density)
    if y is not None:
      interfacial.append([y,z])
    y = fitDensity(ge,density)
    if y is not None:
      interfacial.append([y,z])
    """
    """
  # pyplot.show()
  interfacial = np.array(interfacial)

  return interfacial

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
  # print np.sqrt(pcov[1,1]*(-1/(popt[0]*np.sin(angle)*180/np.pi))**2+pcov[0,0]*(popt[1]/(popt[0]**2*np.sin(angle)*180/np.pi))**2)
  # pyplot.errorbar(interface[:,0], z(interface[:,0], popt[0], popt[1]),yerr=np.sqrt(pcov[0,0]),fmt='yo')
  # pyplot.plot(interface[:,0], z(interface[:,0], popt[0], popt[1]))

  return angle, popt[1], R(*popt).mean()

def main():
  
  # Find the interface for each frame in a trajectory
  trajFile, outFile, nSteps = getArgs(argv)

  t0 = 0
  dt = 10000
  n = 99
  time = dt*n + t0

  inp = None
  try:
    inp = open(trajFile,'r')
    traj = readTraj(inp,['1'],time)
  except IOError:
    if inp is not None:
      inp.close()
    print 'Error reading trajectory file'
  else:
    # read file (need to add a for loop to do every frame)
    frames = []
    traj = combineFrames(traj,nSteps*dt)
    for step in sorted(traj):
      atoms = np.array(traj[step], dtype = float)

      # removes vaporized chains
      # mask = (atoms[:,2] < 20) & (abs(atoms[:,1]) < 80)
      # atoms = atoms[mask]
      # make box
      box = np.array([[np.min(atoms[:,i])-0.01,np.max(atoms[:,i])+0.01] for i in range(3)])
      box[1,:] += [-7,7]

      """
      fig = pyplot.figure()
      ax = fig.add_subplot(1,1,1)
      pyplot.plot(atoms[:,1],atoms[:,2],'o')
      pyplot.ylim(box[2,:])
      pyplot.xlim(box[1,:])
      pyplot.show()
      """
      
      # find interfacial atoms
      interface = findCylinderInterface(atoms, box)
      angle, z0, R = contactAngle(interface[interface[:,1]-box[2,0] > 5], box)
      frames.append([step,angle])
      
    """
      pyplot.plot(interface[:,0],interface[:,1]-box[2,0],'go')
      circ = pyplot.Circle((0,z0), radius=R, fill=False)
      fig.gca().add_artist(circ)
      pyplot.ylabel(r'Height, z ($\sigma$)')
      pyplot.xlabel(r'Width, y ($\sigma$)')
      pyplot.title(r't = %d, $\theta$ = %f' % (step, angle))
      ax.set_aspect('equal','datalim')
      pyplot.show()

    """
    frames = np.array(frames)
    pyplot.figure()
    pyplot.plot(frames[:,0]*0.002,frames[:,1],'bo')
    pyplot.ylabel(r'Contact Angle, $\theta$ (deg)')
    pyplot.xlabel(r'time, t ($\tau_{LJ}$)')
    pyplot.title(r'$\epsilon$ = 0.50')
    # pyplot.ylim((0,80))
    pyplot.show()
  
    print 'Finished analyzing trajectory'

if __name__ == '__main__':
  main()
