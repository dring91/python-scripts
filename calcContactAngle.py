import numpy as np
from sys import argv, exit
from scipy.optimize import curve_fit, leastsq
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from conf_tools import readTrj, readFrame
from lmfit import Model, Parameters

#def combineFrames(file, read_file=readFrame, nCombine=1):
#  
#  time, atoms = readFrame(file)
#  l = len(atoms)
#  combined = np.zeros((l*nCombine,4))
#  combined[0:l] = atoms
#  for n in range(nCombine-1):
#    _, atoms = read_file(file)
#    combined[l*(n+1):l*(n+2)] = atoms
#    
#  return time, combined
# 
#def fitDensity(op, z, density):
#  
#  if op(1,2):
#    guess = [0.8,-40,-1]
#    fstr = 'bo'
#  else:
#    guess = [0.8,40,1]
#    fstr = 'go'
#  mask = op(density[:,0],0)
#  if sum(mask) < 3:
#    return None
#
#  try:
#    params, cov = curve_fit(rho,density[mask][:,0],density[mask][:,1],guess)
#  except RuntimeError:
#    print('z = %.2f slice not fit' % z)
#  else:
#    # Save the "half-density" value to its z-slice
#    if 0.7 < params[0] < 1: # cov[1,1] < 10**3 and not np.isinf(cov[1,1]):      
#      # plotRhoY(z, density,mask,fstr,params)
#      return params[1]
#    
#def findCylinderInterface(atoms, box):
#
#  # define resolution of the z (slices) and y (bins) coordinates
#  # There are ~60k atoms per frame
#  # 60,000 / 15 / 50 = 80 atoms per cell on average
#  # for 5 combined frames, the value is 400
#  zRes = 15
#  yRes = 50
#  # sort atoms according to slice (ragged list)
#  slices, zs, dz = partitionAtoms(atoms, box, zRes)
#
#  interfacial = []
#  rhoz = []
#  # iterate over z slices
#  for index, (slice, z) in enumerate(zip(slices,zs)):
#    # for the nth slice, calculate the radial density
#    density = createHistogram(slice, box, yRes, dz)
#    # fit density to sigmoidal function
#    yl = fitDensity(lt,z,density)
#    if yl is not None:
#      interfacial.append([yl,z])
#    yr = fitDensity(ge,z,density)
#    if yr is not None:
#      interfacial.append([yr,z])
#    # if yl is not None and yr is not None:
#      # rhoz.append([z,density[(yl < density[:,0]) & (density[:,0] < yr)][:,1].mean()])
#    # rhoz.append([z,density[density[:,1] > 0][:,1].mean()])
#  # rhoz = np.array(rhoz)
#  interfacial = np.array(interfacial)
#
#  return interfacial # [interfacial[:,1] == rhoz[:,0][rhoz[:,0] < p[1]]]
#
#def contactAngle(interface,box):
#  # fit interfacial atom coordinates to a circle
#  interface[:,1] -= box[2,0]
#  R = lambda y0, z0: np.sqrt((interface[:,0] - y0)**2 + (interface[:,1] - z0)**2)
#  # res = lambda z0: R(z0) - R(z0).mean()
#  res = lambda args: R(*args) - R(*args).mean()
#  popt, pcov = leastsq(res, [interface[:,0].mean(),-interface[:,1].mean()])
#  # print popt
#  # angle = np.arccos(popt[1]/popt[0])*180.0/np.pi
#  angle = np.arccos(-popt[1]/R(*popt).mean())*180.0/np.pi
#
#  return angle, popt[1], R(*popt).mean()

def rho(z, rho, ze, l): return rho/2*(1-np.tanh(2*(z-ze)/l))
 
def main():
  
  parser = ArgumentParser()
  parser.add_argument('-i','--input')
  parser.add_argument('-o','--output')
  parser.add_argument('-n','--steps',type=int)
  parser.add_argument('-r','--radius',type=float)
  args = parser.parse_args()

  frames = []
  with open(args.input,'r') as inp:
    for n in range(args.steps):
      time, box, atoms = readTrj(inp)

      atoms = atoms[atoms[:,1] == 1][:,2:]
      # Find the center of mass in the xy direction and center droplet
      center = atoms[:,1].mean()
      atoms[:,1] -= center

      # find interfacial atoms
      hist, z_edges, y_edges = np.histogram2d(atoms[:,2], atoms[:,1], bins=(15,60))
      centers = lambda x: np.diff(x) + x[:-1]
      ys, zs = centers(y_edges), centers(z_edges)
      density = hist / np.mean(np.sort(hist.flatten())[::-1][:10])
      print(density.shape)

      fig, (hi,lo) = plt.subplots(nrows=2)
      
      interface_model = Model(rho)
      params = Parameters()
      params.add('rho', value=0.85, vary=True, min=0.8, max=1.0)
      params.add('ze', value=0.0, vary=True, min=35, max=45)
      params.add('l', value=5.0, vary=True, min=1.0, max=6)
      
      l, r = ys<=0, ys>0
      for z,d in zip(zs[:1],density[:1]):
        params['ze'].set(value=40) #np.sqrt(args.radius**2-z**2))
        #left = interface_model.fit(-ys[l],z=d[l],params=params)
        right = interface_model.fit(ys[r],z=d[r],params=params,fit_kws={'ftol':1e-12,'maxfev':1000})
        #hi.plot(-ys[l],left.best_fit)
        hi.plot(ys[r],right.best_fit)
        #hi.plot(-ys[l],d[l])
        hi.plot(ys[r],d[r],label='{:.1f}'.format(z))

      handle = lo.pcolor(ys, zs, density)
      fig.colorbar(handle,orientation='vertical')
      fig.tight_layout()
      plt.show()

    print('Finished analyzing trajectory')

if __name__ == '__main__':
  main()
