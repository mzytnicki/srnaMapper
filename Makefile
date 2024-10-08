CC=gcc
CPPC=g++
ifdef DEBUG
	OFLAGS=-O0 -march=native -g
else ifdef SAN
	OFLAGS=-ggdb -fsanitize=address -fno-omit-frame-pointer -static-libgcc -static-libasan
else ifdef PROF
	OFLAGS=-O3 -pg
else ifdef PERF
	OFLAGS=-O3 -DNDEBUG -march=native -fprofile-use -g
else ifdef VTUNE
	OFLAGS=-O3 -DNDEBUG -march=native -g -shared-libgcc
else
	OFLAGS=-O3 -DNDEBUG -march=native -g -msse2 -mavx
endif
CFLAGS+=-Wall -Wextra $(OFLAGS) -std=c11 -D_GNU_SOURCE
UNAME_S:=$(shell uname -s)
FFLAGS+=-lm -lz -pthread
ifeq ($(UNAME_S),Linux)
	FFLAGS+= -lrt
endif
ifdef ICC
	CC=icc
endif

BWA_O=Libs/bwa/utils.o Libs/bwa/kthread.o Libs/bwa/kstring.o Libs/bwa/ksw.o Libs/bwa/bwt.o Libs/bwa/bntseq.o Libs/bwa/bwa.o Libs/bwa/bwamem.o Libs/bwa/bwamem_pair.o Libs/bwa/bwamem_extra.o Libs/bwa/malloc_wrap.o Libs/bwa/QSufSort.o Libs/bwa/bwt_gen.o Libs/bwa/rope.o Libs/bwa/rle.o Libs/bwa/is.o Libs/bwa/bwtindex.o Libs/bwa/bwashm.o Libs/bwa/bwase.o Libs/bwa/bwaseqio.o Libs/bwa/bwtgap.o Libs/bwa/bwtaln.o Libs/bwa/bamlite.o Libs/bwa/bwape.o Libs/bwa/kopen.o Libs/bwa/pemerge.o Libs/bwa/maxk.o Libs/bwa/bwtsw2_core.o Libs/bwa/bwtsw2_main.o Libs/bwa/bwtsw2_aux.o Libs/bwa/bwt_lite.o Libs/bwa/bwtsw2_chain.o Libs/bwa/fastmap.o Libs/bwa/bwtsw2_pair.o

BWA_C=$(BWA_O:.o=.c)

TARGETS=srnaMapper

all: $(TARGETS)

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

srnaMapper: srnaMapper.c $(BWA_O)
	$(CC) $(CFLAGS) -Wall -Wextra -o srnaMapper.o -c srnaMapper.c
	$(CC) $(CFLAGS) -Wall -Wextra -o srnaMapper srnaMapper.o $(BWA_O) $(FFLAGS)

clean:
	rm -f $(TARGETS) # $(BWA_O)
