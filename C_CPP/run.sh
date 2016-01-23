#!/bin/bash

# ./atom_files.out ../../amits_data/5mers_infiltrate.xyz
# ./msd.out ../../amits_data/5mers_infiltrate.xyz
valgrind --tool=memcheck --leak-check=full ./atom_files.out ../../amits_data/5mers_infiltrate.xyz
