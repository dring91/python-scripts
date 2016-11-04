#!/usr/bin/python

import numpy as np
from sys import argv, exit
from getopt import *
from conf_tools import readFrame, readTrj

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

def findHeight(data, cut_off):
  height = data[0,0]
  prev = data[0]
  for row in data[1:]:
    if row[1] > cut_off >= prev[1]:
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
  hist[:,1] /= 100*100*(1-lower)
  
  return hist

def writeHeader(filename,header):
  # Implicitly assumes a 1 line header
  with open(filename, 'w') as file:
    file.write(header)   
  
def writeHeight(filename,height,mode):
  with open(filename, 'a') as file:
    file.write('%d %f %f %f\n' % tuple(height))

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

  densityFile = path+'density_'+outFile
  heightFile = path+'height_'+outFile
  
  title = '# Infiltration height\n'
  writeHeader(densityFile,title)
  writeHeader(heightFile,title+'# time  cut = 0.85  cut = 0.99  polymer\n')

  mode = 'a'
  nBins = 200

  with open(trajFile, 'r') as file:
    for n in range(nFrames):
      t, box, frame = readTrj(file)

      particles = frame[frame[:,1] == 3][:,4]
      polymers = frame[frame[:,1] == 1][:,4]

      bounds = np.array([particles.min(), particles.max()])
      polymers -= bounds[0]
      polymers /= bounds[1]
      polyLayer = polymers.min()

      density = makeHistogram(polymers,nBins,polyLayer)
      density[:,1] /= bounds[1]
      writeDensity(densityFile,density,t,mode)

      density[:,1] = np.cumsum(density[:,1])/np.sum(density[:,1]*bounds[1])

      height85 = findHeight(density,0.85)
      height99 = findHeight(density,0.99)
      height = [t,height85,height99,polyLayer]
      writeHeight(heightFile,height,mode)
          
  print('Finished analyzing trajectory')

if __name__ == "__main__":  
  main()
