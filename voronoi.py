import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
from conf_tools import readTrj
import argparse
import matplotlib.pyplot as plt

def tetrahedron(point1, point2, point3, point4):

  # find the length of each side of the base
  diff1 = point1-point2
  diff2 = point1-point3
  diff3 = point2-point3
  v_a = diff1.dot(diff1)
  v_b = diff2.dot(diff2)
  v_c = diff3.dot(diff3)
  a = np.sqrt(v_a)
  b = np.sqrt(v_b)
  c = np.sqrt(v_c)

  # calculate the semi-perimeter
  s = (a+b+c)*0.5
  
  # calculate the area with Heron's formula
  heron = np.sqrt(s*(s-a)*(s-b)*(s-c))

  normal = np.cross(v_a,v_b)
  height = (normal/np.sqrt(normal.dot(normal))).dot(point4-point3)

  return heron * height / 3.0

def main():

  """ a script to calculate the density with voronoi diagrams """

  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input")
  parser.add_argument("-s", "--steps")
  args = parser.parse_args()

  with open(args.input, "r") as file:
    for s in (args.steps):
      time, box, atoms = readTrj(file)

  atoms = atoms[np.argsort(atoms[:,0])]
  polymers = atoms[atoms[:,1] == 1][:,2:]
  particles = atoms[atoms[:,1] >= 3][:,2:]
  polymers[:,2] -= particles[:,2].min()
  polymers = polymers[polymers[:,2] >= 0]

  diagram = Voronoi(polymers[:,:2])

  #plt.plot(polymers[:,0], polymers[:,1], '.')
  #plt.plot(diagram.vertices[:,0], diagram.vertices[:,1], '*')
  #plt.xlim(-7,7); plt.ylim(-7,7)
  #plt.show()

  # Select a region and divide it into tetrahedra
  # calculate the region volume with Heron's formula and tetrahedron formula

  voronoi_plot_2d(diagram)
  plt.show()

  #vol = tetrahedron(*points)

if __name__ == '__main__':
  main()
