#!/bin/bash
if [[ $OMPI_COMM_WORLD_RANK == 0 ]]; then
        #nsys profile --trace=cuda,nvtx,mpi,openacc --cuda-memory-usage=true --stats=true $*
	nsys profile --trace=cuda,nvtx,mpi,openacc $*
else
	$*
fi
