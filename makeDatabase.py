import numpy as np
from conf_tools import *
import sqlite3 as sql
import argparse

def main():
  """ script to make a database """

  # parse commandline arguments
  parser = argparse.ArgumentParser()
  parser.add_argument('input')
  args = parser.parse_args()

  with open(args.input,'r') as file:
    box, atoms, bonds = tupleConf(file)

  # Get code to work and then set as context manager
  name = args.input.split('.')[0] + '.db'
  with sql.connect(name) as db:
    field = db.cursor()
    # create and populate atom table
    field.execute('''CREATE TABLE atoms
                     (atom_id integer, mol_id integer, atomtype integer,
                      x real, y real, z real, image_x integer, image_y integer, 
                      image_z integer)''')
    field.executemany('INSERT INTO atoms VALUES (?,?,?,?,?,?,?,?,?)', atoms)
  with sql.connect(name) as db:
    field = db.cursor()
    # create and populate bond table
    field.execute('''CREATE TABLE bonds
                     (bond_id integer, bondtype integer, bond_start integer,
                      bond_end integer)''')
    field.executemany('INSERT INTO bonds VALUES (?,?,?,?)', bonds)
  with sql.connect(name) as db:
    field = db.cursor()
    # create and populate box table
    field.execute('CREATE TABLE box (min real, max real)')
    field.executemany('INSERT INTO box VALUES (?,?)', box)

if __name__ == '__main__':
  main()
