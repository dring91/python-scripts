#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
char **split(char *line, const char *delimiter) {
  char *token = (char *) malloc(sizeof(char *));
  char **tokens = (char **) malloc(sizeof(char **));
  token = strtok(token, delimiter);
  tokens[0] = tokens;
  while (token != NULL) {
    token = strtok(NULL, delimiter);
    (*tokens)++ = tokens;
  }
  return tokens;
}
*/

char *strip(char *line) {
  //check first character for whitespace
  //replace first newline with '\0' and then stop
  char **p = &line;
  while (**p == ' ') {
    (*p)++;
  }
  for (; **p != '\0'; (*p)++) {
    if (**p == '\n')
      **p = '\0';
  }
  return line;
}

//char *readLine(FILE *fp) {
//  size_t bufsize = 10;
//  char buffer[bufsize];
//  // Must include the size of line
//  char 
//
//  while (fgets(buffer, bufsize, fp) != NULL) {
//    if (strchr(buffer,*"\n") == NULL) {
//      line = strcat(line,buffer);
//    } else {
//      line = strcat(line,buffer);
//      break;
//    }
//  }
//  return line;
//}

int 
readAtoms(char  *filename, 
          char   atype, 
           int   nAtoms, 
           int   nDims, 
        double **atoms) {

  size_t size = 50;
  int i=0;

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

void
freeArray(double **array, int size) {
  for (int i=0; i < size; i++)
    free(array[i]);
  free(array);
}

int 
main(int   argc, 
    char **argv) {

  char *filename;
  if (argc == 1) {
    fprintf(stderr,"No cmdline args\n");
    exit(1);
  } else {
    filename = argv[1];
  }
  int nAtoms = 10; // Don't know ahead of time
  int nDims = 3;
  int nLines;
  char atype = '3';

  // allocate atoms array
  double **atoms = (double **) malloc(nAtoms * sizeof(double *)); // [rows][cols]
  if (atoms != NULL) {
    for (int i=0; i < nAtoms; i++)
      atoms[i] = calloc(nDims, sizeof(double));
  }

  // read in atom coordinates and return the number of lines filled
  nLines = readAtoms(filename, atype, nAtoms, nDims, atoms);

  // Should probably resize array using realloc and setting nAtoms to nLines
  // use a temporary pointer and check for success so that memory can be freed
  double **atoms_realloc = (double **) realloc(atoms, nLines * sizeof(*atoms));
  if (atoms_realloc == NULL) {
    freeArray(atoms, nAtoms);
    return 1;
  }
  atoms = atoms_realloc;

  // Loop through array and perform operations
  for (int i=0; i < nLines; i++) {
    printf("%f, %f, %f\n", atoms[i][0], atoms[i][1], atoms[i][2]);
  }

  // free allocated array
  freeArray(atoms, nLines);

  return 0;

}
