import numpy as np
from sys import argv, exit
from scipy.optimize import curve_fit, leastsq
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from conf_tools import readTrj, readFrame
from copy import copy

def rho(z, rhol, ze, l): return rhol/2*(1-np.tanh(2*(z-ze)/l))

def jacobian(z, rhol, ze, l): 

  return np.array([(1-np.tanh(2*(z-ze)/l))/2, 
                  -z/l*rhol/np.cosh(2*(z-ze)/l)**2,
                  rhol*(z-ze)/l**2/np.cosh(2*(z-ze)/l)**2]).T

def contact_angle(interface, weights, box):
  # fit interfacial atom coordinates to a circle by approximating radius
  interface[:,1] -= box[2,0]
  R = lambda y,z: np.sqrt((interface[:,0] - y)**2 + (interface[:,1] - z)**2)
  # calculate objective function
  obj = lambda args: R(*args) - np.average(R(*args), weights=weights/weights.max())
  (y0, z0), _ = leastsq(obj, (interface[:,0].mean(), interface[:,1].mean()))
  # calculate the angle at the surface
  angle = np.arccos(-z0 / R(y0, z0).mean()) * 180/np.pi

  return angle, (y0, z0), R(y0, z0).mean()
 
def main():
  
  parser = ArgumentParser()
  parser.add_argument('-i','--input')
  parser.add_argument('-o','--output')
  parser.add_argument('-n','--steps',type=int)
  parser.add_argument('-r','--radius',type=float)
  parser.add_argument('--plot',action='store_true')
  args = parser.parse_args()

  with open(args.output, 'w') as file: 
    file.write('## Contact angle time series\n# time angle (y0,  z0)\n')
  #with open(args.output + '-interface', 'w') as file:
  #  file.write('## Interface coordinates for contact angle\n# (y,  z)\n')

  frames = []
  interfaces = []
  with open(args.input,'r') as inp:
    for n in range(args.steps):
      time, box, atoms = readTrj(inp)

      atoms = atoms[atoms[:,1] == 1][:,2:]
      # Find the center of mass in the xy direction and center droplet
      center = atoms[:,1].mean()
      atoms[:,1] -= center

      # find interfacial atoms
      hist, z_edges, y_edges = np.histogram2d(atoms[:,2], atoms[:,1], bins=(20,40), range=((0,40),(-50,50)))
      centers = lambda x: np.diff(x)/2 + x[:-1]
      ys, zs = centers(y_edges), centers(z_edges)
      density = hist / np.mean(np.sort(hist.flatten())[::-1][:10])

      if args.plot: fig, (hi,lo) = plt.subplots(nrows=2)
      
      l, r = ys<=0, ys>0
      radius, rhol = 40, 0.8
      guess=(rhol, 1, 0.1)
      interface = []
      params = []
      for z,d in zip(zs,density):
        try: left,  _ = curve_fit(rho, -ys[l]/radius, d[l], p0=guess, method='lm')
        except RuntimeError: left = copy(guess)
        try: right, _ = curve_fit(rho,  ys[r]/radius, d[r], p0=guess, method='lm')
        except RuntimeError: right = copy(guess)
        if left[0] > 0.84 and left[0] < 0.95: 
          if args.plot:
            hi.plot(-ys[l], rho(-ys[l]/radius, *left),color='b')
            hi.scatter(-ys[l],d[l],color='y',label='{:.1f}'.format(z))
            hi.scatter(left[1]*radius,rhol/2,color='r')
          interface.append([-left[1]*radius,z])
          params.append(left)
        if right[0] > 0.84 and right[0] < 0.95: 
          if args.plot:
            hi.plot(ys[r], rho(ys[r]/radius, *right),color='b')
            hi.scatter(ys[r],d[r],color='y',label='{:.1f}'.format(z))
            hi.scatter(right[1]*radius,rhol/2,color='r')
          interface.append([right[1]*radius,z])
          params.append(right)
      params = np.array(params)
      interfaces.append(interface)
      interface = np.array(interface)

      angle, (y0, z0), r = contact_angle(interface, 1/params[:,2], box)

      if args.plot:
        circle = plt.Circle((y0,z0), radius=r, fill=False)
        handle = lo.pcolor(y_edges, z_edges, density)
        cbar = fig.colorbar(handle,orientation='vertical')
        cbar.set_label(r'$\rho\sigma^{3}$')
        lo.scatter(interface[:,0], interface[:,1],color='r')
        lo.add_patch(circle)
        lo.set_aspect('equal')
        hi.set_xlabel(r'$y/\sigma$')
        hi.set_ylabel(r'$\rho\sigma^{3}$')
        lo.set_xlabel(r'$y/\sigma$')
        lo.set_ylabel(r'$z/\sigma$')
        fig.tight_layout()
        plt.show()

      with open(args.output, 'a') as file: 
        file.write('{:d} {:.5f} {:.5f} {:.5f}\n'.format(time, angle, center, z0))
      #with open(args.output + '-interface', 'ab') as file:
      #  np.savetxt(file, interface, fmt='%f %f', comments='', footer=' ')

    np.save(args.output + '-interface.npy', np.array(interfaces))

    print('Finished analyzing trajectory')

if __name__ == '__main__':
  main()
