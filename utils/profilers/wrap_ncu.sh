#!/bin/bash

TAG_TIMESTAMP=`date +"%G%m%d-%H%M%S"`

if [[ $OMPI_COMM_WORLD_RANK == 0 ]]; then
        ncu --set full --kernel-regex ${1} -o report--${1}--${TAG_TIMESTAMP} --target-processes all ${*:2}
else
	$*
fi
