import numpy as np

def PBC(dx,x,interval,cum=False):
  mask = (dx < interval[:,0])
  if cum: mask = np.cumsum(mask,axis=1)
  x = x + 2*interval[:,1]*mask
  mask = (dx >= interval[:,1])
  if cum: mask = np.cumsum(mask,axis=1)
  x = x + 2*interval[:,0]*mask
    
  return x
  
def unwrap(coords, box):
  if len(coords.shape) == 3:
    for i,p in enumerate(coords):
      coords[i,:,:2] = PBC(p[:,:2]-p[0,:2],p[:,:2],box[:2,:])
  elif len(coords.shape) == 2:
    coords[:,:2] = PBC(coords[:,:2]-coords[0,:2],coords[:,:2],box[:2,:])
      
  return coords

def makeBox(dims, padding, *args):
  box = np.zeros((dims,2))
  if len(args) == 0:
    print("No atoms available to make box")
  elif len(args) == 1:
    box[:,0] = args[0].min(0)
    box[:,1] = args[0].max(0)
  elif len(args) >= 2:
    box[:,0] = args[0].min(0)
    box[:,1] = args[1].max(0)

  box[:,0] = box[:,0] - padding * 0.5
  box[:,1] = box[:,1] + padding * 0.5
 
  return box
