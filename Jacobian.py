import numpy as np
import matplotlib.pyplot as plt

def jac(z, rhol, rhov, z0, d):

  sech = lambda x: 1/np.cosh(x)

  J = [0.5 - 0.5 * np.tanh(2*(z-z0)/d),
       0.5 * np.tanh(2*(z-z0)/d) - 0.5,
       (rhol - rhov)*z0/d*sech(2*(z-z0)/d)**2,
       (rhol - rhov)/d**2*(z-z0)*sech(2*(z-z0)/d)**2]
  return J

def main():

  z = np.linspace(0,100,100)
  params = [0.8, 0, 50, 2]
  grads = jac(z, *params)

  for grad in grads:
    plt.plot(z, grad)
  plt.show()

if __name__ == '__main__':
  main()
