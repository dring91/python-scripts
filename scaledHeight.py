#!/usr/bin/python

import numpy as np
from sys import argv, exit
from getopt import *

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
 
def interpolate(x, low, high):
  if x == low[0]:
    return low[1]
  elif x == high[0]:
    return high[1]
  elif x < high[0] and x > low[0]:
    slope = (high[1] - low[1])/(high[0] - low[0])
    return slope*(x - low[0]) + low[1]
  else:
    print 'x out of range'
    exit(0)

def findHeight(data, cut_off=0.85, squared=True):
  height = data[0,0]
  prev = data[0]
  for row in data[1:]:
    if row[1] > cut_off >= prev[1]:
      if squared:
        height = interpolate(cut_off, prev[::-1], row[::-1])**2
      else:
        height = interpolate(cut_off, prev[::-1], row[::-1])
      break
    prev = row
  
  return height
 
def makeHistogram(atoms,nBins,lower):
  hist = np.zeros((nBins,2))
  binSize = (1-lower)/float(nBins)
  hist[:,0] = (np.arange(nBins)+0.5)*binSize+lower
  for atom in atoms:
    if lower < atom <= 1:
      bin = int((atom-lower)/binSize)
      try:
        hist[bin,1] += 1
      except IndexError:
        print "Particle outside box"
  hist[:,1] /= 100*100*binSize
  
  return hist

def writeHeader(filename,header):
  # Implicitly assumes a 1 line header
  with open(filename, 'w') as file:
    file.write(header)   
  
def writeHeight(filename,height,mode):
  with open(filename, 'a') as file:
    file.write('%d' % height[0])
    for h in height[1:]:
      file.write(' %f' % h)
    file.write("\n")

def writeDensity(filename,density,time,mode):
  with open(filename, 'a') as file:
    file.write('#  time: %d \n#   z   density\n' % time)
    [file.write('%f %f\n' % tuple(row)) for row in density]
    file.write('\n')

def getArgs(argv):

  try:
    opts, args = getopt(argv[1:],'t:p:o:n:')
  except GetoptError:
    print 'randomWalk.py -t <xyz file> -p <output path> -o <filename> -n <nFrames>'
    exit(2)
  for opt, arg in opts:
    if opt == '-t':
      inFile = arg
    elif opt == '-p':
      path = arg
    elif opt == '-o':
      outFile = arg
    elif opt == '-n':
      nSteps = arg
 
  return inFile, path, outFile, int(nSteps)
  
def main():

  trajFile, path, outFile, nFrames = getArgs(argv)

  densityFile = path+'density_scaled_'+outFile
  heightFile = path+'height_scaled_'+outFile
  
  title = '# Infiltration height for an undersaturated simulation\n'
  writeHeader(densityFile,title)
  writeHeader(heightFile,title+'# time  cut = 0.85  cut = 0.99\n')

  mode = 'a'
  nBins = 200

  with open(trajFile, 'r') as file:
    for n in range(nFrames):
      t, frame = readFrame(file)

      particles = frame[frame[:,0] == 3][:,3]
      polymers = frame[frame[:,0] == 1][:,3]
      partMin = particles.min()
      # Packing begins at zero
      bounds = np.array([polymers.min()-partMin, 
                         particles.max()-partMin])
      polymers -= partMin
      polymers /= bounds[1]

      density = makeHistogram(polymers,nBins,bounds[0]/bounds[1])
      density[:,1] /= bounds[1]
      writeDensity(densityFile,density,t,mode)

      density[:,1] = np.cumsum(density[:,1])/np.sum(density[:,1])
      height85 = findHeight(density,0.85,False)
      height99 = findHeight(density,0.99,False)
      height = [t,height85,height90,height95,height99,bounds[0]/bounds[1]]
      writeHeight(heightFile,height,mode)
          
  print('Finished analyzing trajectory')

if __name__ == "__main__":  
  main()
