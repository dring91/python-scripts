import numpy as np
import argparse
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def eqns(y, t, a, b, c, d):

  u, h = y

  #return [1 - 6.3*u**2/(5.1*h+30) - 8.8*h*u/(5.1*h+30), u]
  return [1 - a*u**2/(c*h+d) - b*h*u/(c*h+d), u]

def main():
  """ Script to integrate a system of first order differential equations """

  rho = 0.8 # sigma ** -3
  mu = 4.3 # epsilon * tau * sigma ** -3  --> tau = sqrt( sigma ** 2 * m / epsilon )
  r = 5 # sigma
  gamma = 0.39 # epsilon * sigma ** -2
  theta = 0.0*np.pi

  P = 2*gamma*np.cos(theta)/r
  a, b, c, d = 1.225*rho/P, \
               8.0*mu/(P*r**2), \
               rho/P, \
               7.0/6.0*r*rho/P
  c = 0

  y0 = [0, 0]
  t = np.linspace(0,300,100)
  solution = odeint(eqns, y0, t, args=(a,b,c,d))

  #LW = np.sqrt(r*gamma*np.cos(theta)/(2*mu)*t)
  LW = np.sqrt(2/b*t)

  plt.plot(t,solution)
  plt.plot(t,LW)
  plt.xlabel('t')
  plt.ylabel('h,u')
  plt.show()

if __name__ == '__main__':
  main()
