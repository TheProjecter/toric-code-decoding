#!/bin/bash

cd qecsim_noisy_concatenated
make clean
make
nice -n 19 ./qecsim
cd ../qecsim_noisy_toric
make clean
make
nice -n 19 ./qecsim
cd ../qecsim_perfect_concatenated
make clean
make
nice -n 19 ./qecsim
cd ../qecsim_perfect_toric
make clean
make
nice -n 19 ./qecsim
