from sys import argv, exit
import numpy as np
import argparse
import string as st

from conf_tools import readConf

def main():
  
  # get commandline arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("-i","--input")
  parser.add_argument("-o","--output")
  parser.add_argument("sections", nargs=argparse.REMAINDER,
                      help='sections to write from configuration file')
  args = parser.parse_args()

  #sections = {'Atoms','Bonds','Velocities','Masses'}
  sections = {section for section in args.sections}

  # Open input file
  header = None
  with open(args.input, 'r') as file:
    # scan through file
    for line in file:
      items = line.split()
      # scan through sections
      for section in sections:
        if len(items) > 0 and section == items[0]:
          # Assign the current section to the header flag
          header = section
          with open('_'.join([header,args.output]), 'w') as otp:
            otp.write('')
          # print header to show write progress
          print header
      # write section out to file
      if header is not None:
        with open('_'.join([header,args.output]), 'a') as otp:
          otp.write(' '.join([st.strip(line),str(len(items)),'\n']))
  
 
if __name__ == "__main__":
  main()
