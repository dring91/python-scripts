#!/usr/python

import os, string
from sys import argv, exit

# check number of inputs
if (len(argv) < 3):
  print('make-traj.py usage: [input.conf] [otp file]')
  exit(0)

# check for input file
if os.path.exists(argv[1]):

  # open input/output files
  with open(argv[1],'r') as inp:
    with open(argv[2],'w') as otp:
      while (inp):
        L = inp.readline().split()
        # get number of line items and write to otp file
        if ( len(L) > 1 and L[1] == 'atoms'):
          ns = int( L[0] )
          print 'Found', ns, 'atoms'
          line = '%d\n' % ns
          otp.write(line + 'Atoms\n')
  
        if ( len(L) >= 1 and L[0] == 'Atoms'):
          print 'Found atoms section'
          print 'Writing atoms...'
  
          # read in current line
          L = inp.readline()
          for i in range(0,ns):
            L = inp.readline().split()
            line = (L[2]) + ' ' + (L[3]) + ' ' + (L[4]) + ' ' + (L[5]) + '\n'
            otp.write(line)
          break

  print('Finished creating xyz trajectory file from input.conf files')                                        
