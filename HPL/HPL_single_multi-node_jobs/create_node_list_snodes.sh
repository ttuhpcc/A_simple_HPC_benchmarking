#!/bin/bash

#snodes quanah -r -N
sinfo -N -h -o %N -n $(sinfo -hN -p quanah -t idle | grep -v drain | awk '{printf $1","}')
