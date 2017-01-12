import numpy as np
import argparse

def main():
  """ script to calculate/capture the box size for each frame of a trajectory
      and output a configuration that can be used to draw the box
  """

  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input")
  parser.add_argument("-n", "--frames")
  args = parser.parse_args()

  with open(args.input, "r") as file:
    for n in xrange(args.frames):
      time, box, _ = readTrj(file)

