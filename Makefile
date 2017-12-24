CC = gcc
CCFLAGS = -march=native \
	-lm \
	-Wall \
	-O3 \

DBGFLAGS = -march=native \
	-lm \
	-Wall \
	-g \
	-O0 \

all : sbmpart

debug : sbmpart_dbg

sbmpart : sbmpart.c sbm_partition.c sbm_dyn_arr.c sbm_dyn_csr.c sbm_csr.c sbm_file.c sbm_util.c graph.c sbm_stack.c
	$(CC) -o sbmpart sbmpart.c sbm_partition.c sbm_dyn_arr.c sbm_dyn_csr.c sbm_csr.c sbm_file.c sbm_util.c graph.c sbm_stack.c $(CCFLAGS)

sbmpart_dbg : sbmpart.c sbm_partition.c sbm_dyn_arr.c sbm_dyn_csr.c sbm_csr.c sbm_file.c sbm_util.c graph.c sbm_stack.c
	$(CC) -o sbmpart sbmpart.c sbm_partition.c sbm_dyn_arr.c sbm_dyn_csr.c sbm_csr.c sbm_file.c sbm_util.c graph.c sbm_stack.c $(DBGFLAGS)

clean:
	rm sbmpart
