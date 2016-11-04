import numpy as np
from sys import argv, exit
from getopt import *
from memory_profiler import profile
from conf_tools import readFrame, readTrj

def getArgs(argv):
  
  try:
    opts, args = getopt(argv[1:],'i:o:n:')
  except GetoptError:
    print 'randomWalk.py -i <filename> -o <filename> -n <nFrames>'
    exit(2)
  for opt, arg in opts:
    if opt == '-i':
      trajFile = arg
    if opt == '-o':
      outFile = arg
    if opt == '-n':
      nFrames = arg

  return trajFile, outFile, int(nFrames)

def COM(atoms, chainLength):
  size = atoms.shape
  atoms = atoms.reshape((size[0] / chainLength, chainLength, size[1]))
  
  return atoms.mean(1)

def diffusion(atoms_head, atoms_tail, chainLength):
  
  diff_atomic = atoms_head - atoms_tail
  diff_mol = COM(diff_atomic, chainLength)
  dist2_atomic = diff_atomic**2
  dist2_mol = diff_mol**2
  dist2_atomic_var = dist2_atomic.var()
  dist2_mol_var = dist2_mol.var()
  dist2_atomic = dist2_atomic.mean()
  dist2_mol = dist2_mol.mean()

  return dist2_atomic, dist2_mol, dist2_atomic_var, dist2_mol_var
 
# @profile
def main():
  filename, outFile, nFrames = getArgs(argv)

  # Generate multiple MSDs based on the locations of the polymers
  # 1. Exclude polymers in the bulk
  # 2. Exclude polymers NOT in the bulk
  # 3. Look at polymers in the packing exclusively
  # 
  # The point is to find a reasonable time-scale for polymer relaxation
  # Polymers in the packing have z coordinate > 0 and those near the surface are less than
  # 5 sigma from the surface (most negative polymer value)
  MSD = np.zeros((nFrames-1,5))
  #MSD[:,0] = np.arange(1,nFrames)
  # Calculate atomic and molecular MSD
  chainLength = 25
  pos = 0
  for m in range(nFrames-1):
    with open(filename, 'r') as file:
      file.seek(pos)
      time_m, atoms_m = readFrame(file)
      # add a condition to restrict the polymers based on location (make a mask)
      # mask = np.logical_and(atoms_m[:,2] < 0, atoms_m[:,2] > 5*(atoms_m[:,2].min()))
      atoms_m = atoms_m[:,1:]#[mask]
      pos = file.tell()
      for n in range(m+1,nFrames):
        time_n, atoms_n = readFrame(file)
        atoms_n = atoms_n[:,1:]#[mask]
        # handle function calls here
        MSD[n-m-1,1:] += diffusion(atoms_m, atoms_n, chainLength)
        if m == 0:
          MSD[n-m-1,0] = time_n - time_m
  # average MSD array
  MSD[:,1:] /= range(nFrames-1,0,-1)
  # write MSD array to file
  with open(outFile, 'w') as otp:
    otp.write('#  lagtime  MSD\n')
    np.savetxt(otp,MSD,'%0.5f')

if __name__ == '__main__':
  main()
