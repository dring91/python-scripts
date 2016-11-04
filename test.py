import conf_tools

with open('../trajectories/N_25_traj_eps055_1.lammpstrj','r') as file:
  time, box, atoms = conf_tools.readTrj(file)
  print time, box, atoms.shape
