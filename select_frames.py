#!/usr/bin/python

import sys
from conf_tools import *
from pbc_tools import *

MODE = 'r'
X = 0
Y = 1
Z = 2

def main():
  try:
    inFile = sys.argv[1]
  except IndexError:
    print "No Filename given"
    sys.exit()

  nFrames = 20
  nTotal = 200
  nSkip = int(nTotal/nFrames)
  cont = 'yes'

  with open(inFile+'.xyz', MODE) as inp:
    try:
      outFile = sys.argv[2]
    except IndexError:
#       cont = raw_input('Warning: output filename is the same as the input file. This may overwrite data. Continue? [y]/n: ')
#       if cont == 'no' or cont == 'n':
#         sys.exit()
#       else:
      outFile = inFile

    suffix = '_sparse'
    for n in range(nFrames):
      for k in range(nSkip):
        time, frame = readFrame(inp)
      write_xyz(outFile+suffix, time, frame, 'a')

if __name__ == '__main__':
  main()
