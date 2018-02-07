import numpy as np
from argparse import ArgumentParser

from conf_tools import readConf

def main():
  """ Calculate:
        1. whether a polymer is adsorbed
        2. what is its degree of adsorption
        3. size of adsorbed layer
        4. correlation of adsorption with height
  """

  parser = ArgumentParser()
  parser.add_argument('-i','--input')
  parser.add_argument('-o','--output')
  parser.add_argument('-n','--frames',type=int)
  parser.add_argument('-m','--nMon',type=int)
  parser.add_argument('-r','--radius',type=float)
  parser.add_argument('--polymer',type=int)
  parser.add_argument('--cylinder',type=int)
  args = parser.parse_args()

  with open(args.input, 'r') as file:
    for n in range(args.frames):
      time, box, atoms = readConf(file)    

      atoms = atoms[atoms[:,0].argsort()]
      polymers = atoms[atoms[:,1] == args.polymer][:,2:]
      cylinder = atoms[atoms[:,1] == args.cylinder][:,2:]

      polymers = polymers[polymers[:,2] > cylinder[:,2].min()]
      polymers[:,:2] -= cylinder[:,:2].mean(0)
      if len(polymers) == 0: continue
      polymers = polymers.reshape((-1,args.nMon,3))

      layer = 1.5
      adsorbed = polymers[:,:,0]**2 + polymers[:,:,1]**2 > (args.radius - layer)**2
      chains = adsorbed.sum(1) 
          

if __name__ == '__main__':
  main()
