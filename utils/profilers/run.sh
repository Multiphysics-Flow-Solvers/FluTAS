rm -rf data/* out && mpirun -np 1 --map-by socket:PE=5 --mca btl ^openib --report-bindings ./wrap_nsys.sh ../../../src/flutas 2>&1 | tee out 
