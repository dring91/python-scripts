import numpy as np
from sys import argv, exit
from conf_tools import *
import matplotlib.pyplot as plt

def main():
  try:
    logfile = argv[1]
    num = argv[2]
    nFrames = int(argv[3])
  except IndexError:
    print 'missing input'
    exit(1)

  pressure = np.loadtxt(logfile, float, usecols=(1,2,3,6), ndmin=2)
  nSkip = 0
  nBins = 100
  dBin = pressure[0,3] / float(nBins)
  # nFrames = 100
  Lz = pressure[0,3]

  ave = np.zeros(nFrames)
  std = np.zeros(nFrames)
  for i in range(nFrames):
    ave[i] = 0.5*Lz*(pressure[:(i+1),2].mean()-0.5*(pressure[:(i+1),1].mean()+pressure[:(i+1),0].mean()))
    std[i] = np.sqrt(pressure[:(i+1),0].var()*(-0.25*Lz)**2 + 
                     pressure[:(i+1),1].var()*(-0.25*Lz)**2 +
                     pressure[:(i+1),2].var()*(0.5*Lz)**2)/np.sqrt(i)

  # plt.errorbar(range(1,nFrames+1),ave,yerr=std/2)
  # plt.show()

  with open("surf-ten_N%s" % num,"w") as otp:
    for line in zip(ave,std):
      otp.write("%f %f\n" % line)

  gamma = 0.5*Lz*(pressure[:,2]-0.5*(pressure[:,0]+pressure[:,1]))
  print gamma[nSkip:].mean(), gamma[nSkip:].std()/np.sqrt(nFrames)

if __name__ == '__main__':
  main()
