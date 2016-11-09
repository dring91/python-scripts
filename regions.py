import numpy as np

def inRange(point, interval=[0,1], left=True, right=False):
  
  if point > interval[0] and point < interval[1]:
    return True
  elif left and point == interval[0]:
    return True
  elif right and point == interval[1]:
    return True
  else:
    return False
      
def inBox(point, box=[[0,1],[0,1],[0,1]], left=True, right=False):
  
  for i in range(len(point)):
    if inRange(point[i],box[i]) == False:
      return False
  
  return True

def inSphere(point, R=[0,1], inner=True, outer=True, hemisphere=False):
  
  r = np.sqrt(np.dot(point,point))
  if hemisphere and inRange(r,R,inner,outer) and point[2] > 0:
    return True
  elif inRange(r,R,inner,outer):
    return True
  
  return False

def inCylinder(point, R, L, origin=[0,0,0], outer=True, inner=True):

  R = [0.0,R]
  L = [-0.5*L,0.5*L]
  point = np.array(point)
  origin = np.array(origin)
  r = np.sqrt(np.dot(point[1:]-origin[1:],point[1:]-origin[1:]))
  if inRange(r,R,inner,outer) and inRange(point[0]-origin[0],L):
    return True

  return False
