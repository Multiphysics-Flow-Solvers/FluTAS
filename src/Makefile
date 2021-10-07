# Makefile
ifndef ARCH
ARCH := generic
endif

ifndef APP
APP := basic
endif

include targets/target.$(ARCH)
include apps/app.$(APP)

POST_PATH = apps/postp.$(APP) -I./postprocessing

.PHONY: clean

.SUFFIXES :
.SUFFIXES : .o .c .f90

.f90.o:
	$(FC) $(FFLAGS) $(DFLAGS) $(CPP_FLAGS) -c -I$(POST_PATH) -o $@ $<

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

OBJS = bound.o \
chkdiv.o \
chkdt.o \
common_mpi.o \
correc.o debug.o \
fft.o \
fftw.o \
fillps.o \
initflow.o \
initgrid.o \
initmpi.o \
initsolver.o \
load.o \
output.o \
param.o \
rk.o \
source.o \
sanity.o \
2decomp/decomp_2d.o \
2decomp/io.o \
types.o

###########################################################
######## TO CUSTOMIZE BASED ON APPLICATION PHYSICS ########
###########################################################

ifeq ($(TURB_FORCING),1) 
	DFLAGS += -D_TURB_FORCING
endif

ifeq ($(DO_POSTPROC),1) 
	OBJS += postprocessing/post.o
# OBJS += postprocessing/spectra.o
# OBJS += postprocessing/fftw_spectra.o 
# OBJS += postprocessing/post.o 
# OBJS += postprocessing/tagging.o
	DFLAGS += -D_DO_POSTPROC
endif

ifeq ($(INIT_MONTECARLO),1) 
	DFLAGS += -D_INIT_MONTECARLO
	OBJS += monte_sphere.o
endif

ifeq ($(CONSTANT_COEFFS_POISSON),1) 
	DFLAGS += -D_CONSTANT_COEFFS_POISSON 
	OBJS += solver.o solver_cpu.o
else 
	OBJS += solver_vc.o
	LDFLAGS += -lHYPRE
endif

ifeq ($(HEAT_TRANSFER),1) 
	DFLAGS += -D_HEAT_TRANSFER -D_MOM_QUICK
ifeq ($(BOUSSINESQ),1)
        DFLAGS += -D_BOUSSINESQ
endif
	OBJS += gradls.o moms.o rks.o cmpt_divth.o
endif

ifeq ($(USE_VOF),1) 
	DFLAGS += -D_USE_VOF
	OBJS += vof.o mom.o funcs.o
else
	OBJS += mom.o moms.o
endif

ifeq ($(USE_NVTX),1)  
	DFLAGS += -D_USE_NVTX
	OBJS += nvtx.o
  LDFLAGS += $(NVTX_LIB)
endif

###########################################################

TARGET = flutas

all: $(TARGET)

main.o : apps/$(MAIN_SRC) $(OBJS)
	$(FC) $(FFLAGS) $(DFLAGS) $(CPP_FLAGS) -I$(POST_PATH) -c -o main.o $<

$(TARGET): $(OBJS) main.o
	$(FC) $(FFLAGS) $(DFLAGS) $(OBJS) $(LDFLAGS) main.o -o $(TARGET)

clean:
	rm -rf *.o *.mod *dSYM $(TARGET)

include make.deps
