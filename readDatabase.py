import numpy as np
import sqlite3 as sql
import argparse

def main():
  """ sample script to read a database """

  # parse commandline arguments
  parser = argparse.ArgumentParser()
  parser.add_argument('input')
  args = parser.parse_args()

  with sql.connect(args.input) as db:
    field = db.cursor()
    # create and populate atom table
    field.execute('''SELECT * FROM atoms WHERE atomtype=?''', ('1',))
    polymers = np.asarray(field.fetchall())
    field.execute('''SELECT * FROM box''')
    box = np.asarray(field.fetchall())

if __name__ == '__main__':
  main()
