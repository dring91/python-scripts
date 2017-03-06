from sys import argv, exit
import numpy as np
import argparse
import matplotlib.pyplot as plt

from conf_tools import readConf, write_conf, write_xyz, write_traj
from pbc_tools import unwrap, PBC

def main():
  
  # get commandline arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("-i","--input")
  parser.add_argument("-o","--output")
  parser.add_argument("-d","--descrip",help="Header for file describing what it contains")
  parser.add_argument("-l","--length",help="chain length of polymer molecule",type=int)
  parser.add_argument("-t","--atomtype",help="atomtype ID for layer to unwrap",type=int)
  parser.add_argument("--no-image-flags",dest="entries",action="store_const",const=6,
                      default=9,help="switch to adjust number of entries for data flags")
  args = parser.parse_args()

  with open(args.input, 'r') as file:
    for line in file:
      if line.split()[0] == 'Atoms':
        
  
 
if __name__ == "__main__":
  main()
