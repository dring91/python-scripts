#!/usr/bin/python

import numpy as np

def PBC(dx,x,interval):
  # isn't this redundant?
  #   x = x * mask + x * np.logical_not(mask)
  mask = (dx < interval[:,0])
  x = x*mask + 2*interval[:,1]*mask + x*np.logical_not(mask)
  mask = (dx >= interval[:,1])
  x = x*mask + 2*interval[:,0]*mask + x*np.logical_not(mask)
    
  return x
  
def unwrap(coords, box):
  for i,p in enumerate(coords):
    coords[i,:,:2] = PBC(p[:,:2]-p[0,:2],p[:,:2],box[:2,:])
      
  return coords

def makeBox(dims, *args):
  box = np.zeros((dims,2))
  if len(args) == 0:
    print "No atoms available to make box"
  elif len(args) == 1:
    box[:,0] = args[0].min(0)
    box[:,1] = args[0].max(0)
  elif len(args) >= 2:
    box[:,0] = args[0].min(0)
    box[:,1] = args[1].max(0)
 
  return box
