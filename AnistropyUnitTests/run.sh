#!/bin/bash
#print_args.sh
N="$@"
echo "Number of points:" $N
make clean all
./test_anisotropy $N

