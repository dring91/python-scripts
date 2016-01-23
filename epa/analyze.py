#!/usr/bin/python

import numpy as np
from sys import argv
from sys import stdout
from importlib import import_module

def loadPlugIn(moduleName):
  # should load the module containing the specific calculations needed
  # import moduleName as operation
  operation = import_module(moduleName)
  return operation

def readFile(file):
  line = None
  data = []
  while ( line != '' ):
    line = file.readline()
    L = line.split()
    if len(L) > 0 and L[0] != '#':
      L = np.array(L)
      data.append(L)

  return data

def writeNums(file, nums, whitespace='\n'):
  try:
    for num in nums:
      writeNums(file, num, whitespace)
  except TypeError:
    file.write(str(nums))
    file.write(' ')
  else:
    file.write(whitespace)

def writeData(file, nums, header='# col1 col2 ...\n'):
  # Call polymorphic write function and catch when it runs out of dimensions 
  file.write(header)
  writeNums(file, nums, '\n')

def printData(nums):
  writeData(stdout, nums, '')
      
def main():
  
  plug_in = loadPlugIn(argv[2])
  print '# Plug-in loaded'
  
  with open(argv[1], 'r') as file:
    data = readFile(file)
  data = np.array(data)
  print '# data loaded'

  try:
    args = argv[4:]
  except IndexError:
    args = []
  result = plug_in.command(data,*args)
  print '# computation complete'

  header = '# test\n'
  # May want to ask in future if the user wants to overwrite the file
  mode = 'w'

  try:
    outFile = argv[3]
  except IndexError:
    printData(result)
  else:
    with open(outFile, mode) as file:
      writeData(file, result, header)
  print '# data written'

if __name__ == "__main__":
  main()
