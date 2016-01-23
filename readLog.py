#!/bin/python

from sys import argv
import numpy as np

def readFile(filename, nFrames):
  with open(filename, 'r') as file:
    i = 0
    for line in file:
      L = line.split()
      if len(L) > 0 and L[0] == "Step":
        data = np.zeros((nFrames,len(L)))
        labels = L
      try:
        data[i] = np.array(L,dtype=float)
      except (NameError, ValueError, IndexError):
        pass
      else:
        i += 1

  return labels, data

def writeFile(filename, mode, labels, data):
  with open(filename, mode) as file:
    file.write("# ")
    [file.write("%s " % l) for l in labels]
    file.write("\n")
    try:
      for row in data:
        [file.write("%f " % col) for col in row]
        file.write("\n")
    except TypeError:
      [file.write("%f " % d) for d in data]
      file.write("\n")
    file.write("\n")

def main():

  # Get the path (and possibly the log name)
  path = argv[1]
  # assume that inp name is log.lammps and otp name is thermo
  inFile = "log.lammps"
  outFile = "thermo"
  nFrames = 201
  # read file by searching for line beginning with "Step"
  # read in the header as well so that the data can be referenced
  labels, data = readFile(path+inFile, nFrames)
  writeFile(path+outFile, "w", labels, data)
  # Calculate and write statistics
  stats = data.mean(0)
  writeFile(path+outFile, "a", labels, stats)

if __name__ == "__main__":
  main()
