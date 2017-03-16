-- Open a connection to the database
.open test.db

-- Make one table referencing each frame
create table frames(
  frameid        INTEGER PRIMARY KEY,
  timestep       INTEGER UNIQUE,
  atoms_table_id TEXT UNIQUE,
  box_table_id   TEXT UNIQUE
);

-- Make two tables for each frame, one for the atoms and one for the boxes
create table atoms_0(
  atomid  INTEGER PRIMARY KEY,
  molid   INTEGER,
  x_coord REAL,
  y_coord REAL,
  z_coord REAL
);

create table box_0(
  xlo REAL,
  xhi REAL,
  ylo REAL,
  yhi REAL,
  zlo REAL,
  zhi REAL
);

-- Make link between main table and frame tables
insert into frames(frameid, timestep, atoms_table_id, box_table_id) values(0, 0, 'atoms_0', 'box_0');
