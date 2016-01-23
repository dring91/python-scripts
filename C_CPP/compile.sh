#!/bin/bash

name=$1
gcc -std=c11 -g -Wall ${name}.c -o ${name}.out
