import numpy as np
import matplotlib.pyplot as plt

def rotations(v, angle, axis):

  axes = np.arange(3)
  [i,j] = axes[axes != axis]

  R = np.zeros((3,3))
  R[axis,axis] = 1
  R[i,i] = np.cos(angle) 
  R[j,j] = np.cos(angle) 
  R[i,j] = np.sin(angle) 
  R[j,i] = -np.sin(angle) 

  return R.dot(v)

def ellipsoid(coords, phi=0, theta=0, axes=[1,1,1]):

  coords -= coords.mean(0)
  A = np.eye(3)*np.array(axes,dtype=float)**(-2)

  Ay = rotations(A, theta, 1)
  Ayz = rotations(Ay, phi, 2)

  Ayz = np.einsum('ij,ki->kj',Ayz,coords)
  Ayz = np.einsum('ij,ij',Ayz,coords)

  Ayz -= 1
 
  return Ayz

def main():

  v = np.ones(3)
  theta = np.pi / 2
  phi = np.pi / 4

  u = rotations(v,theta,0)

  w = rotations(u,phi,2)

  arrows = np.zeros((3,4))
  arrows[:,1] = v
  arrows[:,2] = u
  arrows[:,3] = w

  arrows = arrows*np.logical_not(np.isclose(arrows,np.zeros_like(arrows)))
  print arrows

  #plt.plot(arrows[0,:2],arrows[1,:2])
  l = 3
  plt.plot(arrows[1,0:l:(l-1)],arrows[2,0:l:(l-1)])
  plt.show()

if __name__ == '__main__':
  main()
