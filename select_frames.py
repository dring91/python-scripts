from sys import argv,exit
from conf_tools import *
from pbc_tools import *
import argparse

MODE = 'r'

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input", required=True)
  parser.add_argument("-o", "--output", required=True)
  parser.add_argument("-t", "--step")
  args = parser.parse_args()

  suffix = '_sparse'
  nFrames = 200

  with open(args.input+'.xyz', MODE) as inp:
    for n in range(nFrames):
      time, frame = readFrame(inp)
      if time % args.step == 0:
        write_xyz(args.output+suffix, frame, time, 'a')

if __name__ == '__main__':
  main()
