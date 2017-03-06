import numpy as np

class Trajectory(object):

  def __init__(self,file):
    """ initialize class """

class Frame:

  def __init__(self):
    """ initialize a trajectory frame """
    self.nMolecules
    self.nAtoms
    self.nBonds
    self.nTypes
    self.time
         
    self.box = np.zeros((3,2))
    self.atoms = np.zeros((nAtoms, 6))
    self.bonds = np.zeros((nBonds, 4))

  def load(self,file):
    """ load frame into class """
