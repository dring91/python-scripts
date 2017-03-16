import numpy as np
from sys import exit
import argparse
from conf_tools import *
import sqlite3 as sql

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-t", "--text", help='text file to convert to database')
  parser.add_argument("-d","--database")
  args = parser.parse_args()

  # connect to the sqlite database
  connection = sql.connect(args.database,isolation_level="DEFERRED")

  # create main table
  try: connection.execute('''create table frames(
                               frameid        INTEGER PRIMARY KEY,
                               timestep       INTEGER UNIQUE,
                               num_atoms      INTEGER,
                               xlo            REAL,
                               xhi            REAL,
                               ylo            REAL,
                               yhi            REAL,
                               zlo            REAL,
                               zhi            REAL,
                               atoms_table_id TEXT UNIQUE
                             );''')
  except sql.DatabaseError: print 'table "frames" already created'
  finally: connection.commit()

  with open(args.text, "r") as file:
    for i, (time, n, box, atoms) in enumerate(iTrj(file)):

      print i

      # create frame tables
      try: connection.execute('''create table atoms_%d(
                                   atomid  INTEGER PRIMARY KEY,
                                   typeid  INTEGER,
                                   x_coord REAL,
                                   y_coord REAL,
                                   z_coord REAL
                                 );''' % i)
      except sql.DatabaseError: print 'table "atoms_%d" already created' % i
      finally: connection.commit()

      # create table links
      atom_id = 'atoms_%d' % i
      box = box.reshape(6).tolist()
      frame_info = [i, time, n] + box + [atom_id]
      try: connection.execute('insert into frames values(?,?,?,?,?,?,?,?,?,?);', 
                              tuple(frame_info)
                             )
      except sql.IntegrityError: print 'ignoring duplicate primary key in "frames"'
      finally: connection.commit()

      # Insert atoms into atom table
      insert_str = 'insert into atoms_%d values(?,?,?,?,?)' % i
      try: connection.executemany(insert_str,atoms)
      except sql.IntegrityError: print 'ignoring duplicate atom_ids in atoms_%d' % i
      finally: connection.commit()

  connection.close()

if __name__ == '__main__':
  main()
