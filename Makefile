# Build all flutas.* binaries at once

.PHONY: all clean

.NOTPARALLEL:

ARCH ?= generic
USE_FAST_KERNELS ?= 0
USE_NVTX ?= 0

APP_LIST=basic two_phase_ht two_phase_inc_isot two_phase_inc_isot_turb

all: 
	@for idapp in $(APP_LIST); \
	do \
		make -C ./src clean-obj; make -C ./src ARCH=$(ARCH) APP=$${idapp} USE_FAST_KERNELS=$(USE_FAST_KERNELS) USE_NVTX=$(USE_NVTX) -j4 flutas.$${idapp}; \
	done

clean-obj:
	make -C ./src clean-obj

clean: clean
	@for idapp in $(APP_LIST); \
	do \
		make -C ./src ARCH=$(ARCH) APP=$${idapp} clean; \
	done