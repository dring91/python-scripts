#include <stdio.h>
#include <stdlib.h>

// reads the atom coordinates into a 3D array and filters out NP/surface atoms
// array[nFrames][nAtoms][nDims]
int 
readArray(char *filename, 
          char  atype,
           int  nPoly, 
           int  nDims, 
         float *atoms[nPoly][nDims]) {

  size_t size = 50;
  int i=0;
  int nAtoms = 272844;
  nAtoms += nPoly;

  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    printf("File cannot be opened\n");
  } else {
    // Read line-by-line with fgets
    char line[size];
    char words[4][size];
    while (nAtoms-- > 0) {
      if (fgets(line, size, fp) != NULL) {
        // filter lines read and split into variables
        sscanf(line, "%s %s %s %s", words[0], words[1], words[2], words[3]);
        if (*words[0] == atype) {
          for (int j=0; j < nDims; j++) {
            atoms[i][j] = atof( words[j+1] );
          }
          i++;
        }
      }
    }
    fclose(fp);
  }

  return i;
}

//   unwraps the atom coordinates and keeps track of
// how many times an atom has crossed the boundary
void unwrap() {

}

// writes ths final MSD to a file
void writeData() {

}

// calculates the Mean Squared Displacement of the atoms
void MSD() {

}

int main(int argc, char **argv) {

  // read in argv and assign
  char *filename = argv[1];

  // declare variables
  //int nFrames = 200;
  char atype = '1';
  int nAtoms = 350000;
  int nDims = 3;
  float atoms[nAtoms][nDims];

  // read in the data
  readArray(filename, atype, nAtoms, nDims, &atoms);
  for (int i=0; i < 100; i++) {
    printf("%f, %f, %f\n", atoms[i][0], atoms[i][1], atoms[i][2]);
  }

  // unwrap the coordinates at every timestep

  // calculate the MSD

  // write the MSD to a file

  return 0;
}
