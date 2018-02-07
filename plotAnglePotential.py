import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def main():

  fig, ax = plt.subplots()
  mpl.rc('font',size=12)

  theta = np.linspace(0,np.pi,30)

  ax.plot(theta, 5*(theta-2*np.pi/3)**2, label='harmonic')
  ax.plot(theta, 3*(1+np.cos(theta)), label=r'cosine, $K = 3$')

  ax.legend()
  ax.set_xlim(0,np.pi)
  ax.set_xlabel(r'$\theta$',size=12)
  ax.set_ylabel(r'$U(\theta)$',size=12)
  plt.tight_layout()
  plt.show()

if __name__ == '__main__':
  main()
