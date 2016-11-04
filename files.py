# This file defines the xyz class. It is derived from a class called File and is used to make switching between file formats easier.
import numpy as np

class File:

  def __init__(self, filename):
    self.name = filename
    self.nAtoms = 0
    self.time = []
    self.atoms = np.array([])

  # provides context management functionality
  def __enter__(self):
    return self

  def __exit__(self, exc_type, exc_val, traceback):
    file.close()

class XYZ(File):

  def readFile(self):

  def writeFile(self):

class TRJ(File):

  def readFile(self):

  def writeFile(self):
    file.write('ITEM: TIMESTEP\n%d\n' % self.time)
    file.write('ITEM: NUMBER OF ATOMS\n%d\n' % len(frame))
    file.write('ITEM: BOX BOUNDS pp pp ss\n')
    np.savetxt(file, box, '%f %f')
    file.write('ITEM: ATOMS id type xs ys zs\n')
    np.savetxt(file, frame, '%d %d %f %f %f\n')
