#!/bin/bash

cd qecsim_noisy_concatenated
gprof ./qecsim | gprof2dot.py | dot -Tsvg -o noisy_concatenated_cg.svg
cd ../qecsim_noisy_toric
gprof ./qecsim | gprof2dot.py | dot -Tsvg -o noisy_toric_cg.svg
cd ../qecsim_perfect_concatenated
gprof ./qecsim | gprof2dot.py | dot -Tsvg -o perfect_concatenated_cg.svg
cd ../qecsim_perfect_toric
gprof ./qecsim | gprof2dot.py | dot -Tsvg -o perfect_toric_cg.svg
