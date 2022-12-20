#!/bin/bash
#print_args.sh
N="$@"
echo "Number of points:" $N
make clean all OPT=yes OMP=yes
./test_anisotropy $N

