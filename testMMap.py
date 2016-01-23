#!/usr/bin/python

import os
import mmap

def memory_map(filename, access=mmap.ACCESS_WRITE):
  size = os.path.getsize(filename)
  filedesc = os.open(filename, os.O_RDRW)
  return mmap.mmap(filedesc, size, access=access)

def main():

  filename = "../685558.ocoee.nics.utk.edu/N_1_traj_eps100_T50_1.xyz"
  
  # with memory_map(filename) as mm:
  #   index = mm.find("Timestep: ")
  #   print index #, mm[index]

  with open(filename, 'r+b') as f:
    mm = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
    print mm.readline()
