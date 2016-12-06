from sys import argv, exit
import numpy as np
import argparse

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument("-o","--output")
  parser.add_argument("-c","--nChains", type=int)
  parser.add_argument("-m","--nMons", type=int)
  parser.add_argument("-ar", type=float)
  args = parser.parse_args()

  # total number of lj beads in melt
  total = args.nMon * args.nChains

  # approximate melt density and volume
  density = 0.8
  volume = total / density
  
  # calculations for finding box dimensions
  # AR = L / h
  # V = L * L * h
  # solve for thickness: h = cube_root( V / AR^2 )
  # L = h * AR
  #
  # eg: AR = 1 N = 10^5
  #     V = 125,000
  #     h = 50
  #     L = 50
  # eg: AR = 1/3 N = 10^5
  #     V = 125,000
  #     h = 24
  #     L = 72
  # 
  # V = (L / d) * (L / d) * (h / d) * d^3
  # V / d^3 = N = XS * YS * ZS
  # d = cube_root ( N / V )
  #
  # eg: AR = 1 N = 10^5
  #     d = 0.928
  #     XS = YS = L / d = h * AR / d = 54
  #     ZS = N / 54 ^ 2 = 42.87

  # make box
  box = np.zeros((3,2))
  h = (volume / ar**2)**(1/3)
  box[2,1] = h
  box[:2,1] = h * ar
  box -= 0.5 * (box[:,1] + box[:,0])

  # coordinate placement
  # [ X1, Y1, Z1 ]
  # [ X2, Y1, Z1 ]
  # [ X1, Y2, Z1 ]
  # [ X2, Y2, Z1 ]
  # [ X1, Y1, Z2 ]
  # [ X2, Y1, Z2 ]
  # [ X1, Y2, Z2 ]
  # [ X2, Y2, Z2 ]

  # generate particle positions inside box
  # generate too many positions and then truncate the coordinate array
  XS, YS, ZS = 54, 54, 43
  dimensions = 54, 54, 43
  coordinates = np.zeros((total, 3))
  coordinates = np.indices(dimensions)
  coordinates = coordinates.reshape((3,-1)).T
  
if __name__ == '__main__':
  main()
