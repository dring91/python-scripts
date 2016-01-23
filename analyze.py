#!/usr/bin/python

import numpy as np
from sys import argv
from importlib import import_module

# get file name and open file
#  -- analyze file mean output --
#     ^
#     analyze calls the wrapping program
#             file
#             ^
#             file is the name of the data file
#                  ^
#                  mean is an example of an operation script that the wrapper can use
#                       output
#                       ^
#                       output is the location to output to
#
# read file data line by line
#  1. Use the # as a comment
#  2. split all lines and convert to floats
# 
# after each line read, perform an action
# There are 2 main types of actions: 
#  1. line operations
#  2. accumulators (operate on a line, but keep track of a value in between reading
#
# when reading has finished perform final actions and exit
# Final actions may include:
#  1. writing to a file
#  2. manipulating an accumulator

def loadModule(moduleName):
  # should load the module containing the specific calculations needed
  # import moduleName as operation
  operation = import_module(moduleName)
  return operation

def readLine(file):
  # return line
  line = file.readline()
  line = line.split()
  if len(line) > 0 and line[0] != '#':
    line = np.array(line)
  elif len(line) > 0 and line[0] == '#':
    line = None

  if line == None:
    line = readLine(file)
  
  return line

def manipulate(line, accumulator=0):
  # should return whatever operation returns
  pass

def writeData(filename, data, mode, header='# col1 col2 ...\n'):
  with open(filename, mode) as file:
    file.write(header)
    try:
      (nRows, nCols) = data.shape
    except ValueError:
      (nCols,) = data.shape
      for j in range(nCols):
        file.write(str(data[j]))
        file.write(' ')
      file.write('\n')
    else:
      for i in range(nRows):
        for j in range(nCols):
          file.write(str(data[i,j]))
          file.write(' ')
        file.write('\n')

def main():
  
  operation = loadModule(argv[2])
  
  line = None
  data = []
  with open(argv[1], 'r') as file:
    while ( line != [] ):
      line = readLine(file)
      if line != []:
        data.append(line)

  data = np.array(data, dtype=float)
  result = operation.command(data, 'rows')
  result = np.array(result)

  try:
    outFile = argv[3]
  except IndexError:
    print 'No output filename given'

  # header = '# unweighted  weighted\n'+
  # header = '# insertion method\n'
  # header = '# polymer chain method\n'
  writeData(argv[3], result, 'a', header)

if __name__ == "__main__":
  main()
