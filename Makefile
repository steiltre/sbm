CC = gcc
CCFLAGS = -march=native \
	-lm \
	-Wall \
	-O3 \
	-fopenmp

DBGFLAGS = -march=native \
	-lm \
	-Wall \
	-g \
	-O0 \
	-fopenmp

LDLIBS = -lpcg_random

all : sbmpart

debug : sbmpart_dbg

sbmpart : sbmpart.c sbm_partition.c sbm_dyn_arr.c sbm_dyn_csr.c sbm_csr.c sbm_file.c sbm_util.c graph.c sbm_stack.c
	$(CC) -o sbmpart sbmpart.c sbm_partition.c sbm_dyn_arr.c sbm_dyn_csr.c sbm_csr.c sbm_file.c sbm_util.c graph.c sbm_stack.c $(CCFLAGS) $(LDLIBS)

sbmpart_dbg : sbmpart.c sbm_partition.c sbm_dyn_arr.c sbm_dyn_csr.c sbm_csr.c sbm_file.c sbm_util.c graph.c sbm_stack.c
	$(CC) -o sbmpart sbmpart.c sbm_partition.c sbm_dyn_arr.c sbm_dyn_csr.c sbm_csr.c sbm_file.c sbm_util.c graph.c sbm_stack.c $(DBGFLAGS) $(LDLIBS)

clean:
	rm sbmpart
