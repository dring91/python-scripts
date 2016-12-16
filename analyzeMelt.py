import numpy as np
from conf_tools import *
from pbc_tools import *
import argparse


def main():
  
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input")
  parser.add_argument("-o", "--output")
  parser.add_argument("-n", "--frames", type=int)
  args = parser.parse_args()

  with open(args.output, "w") as out:
    out.write("time bond_length variance end-to-end variance\n")

  with open(args.input, "r") as file:
    for n in range(args.frames):
      time, box, frame = readTrj(file)

      frame = frame[np.argsort(frame[:,0])][:,2:]
      frame = frame.reshape((400,100,3))
      frame = unwrap(frame,box)

      bond_length = 0.97
      overlap = 0
      for c in range(1,100):
        diff = frame[:,c]-np.swapaxes(frame[:,:c],0,1)
        overlap += (np.sqrt(np.einsum('ijk,ijk->ji',diff,diff)) < bond_length).sum()
      print overlap

      diff = frame[:,1:] - frame[:,:-1]
      bonds = np.sqrt(np.einsum('...i,...i',diff,diff))

      R = diff.sum(1)
      R = np.sqrt(np.einsum('...i,...i',R,R))

      with open(args.output, "a") as out:
        out.write("%d %f %f %f %f\n" % (time, bonds.mean(), bonds.var(), R.mean(), R.var()))

if __name__ == '__main__':
  main()
