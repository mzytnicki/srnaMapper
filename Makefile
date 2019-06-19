ifdef DEBUG
	OFLAGS=-O3 -DNDEBUG -march=native -g
else ifdef PROF
	OFLAGS=-O3 -pg
else
	OFLAGS=-O3 -DNDEBUG -march=native -g
	#OFLAGS=-O3 -DNDEBUG -march=native -fprofile-use -g
endif
CPPC=clang++
CC=gcc
CFLAGS=-Wall -Wextra $(OFLAGS) -pthread
LFLAGS=-lrt

BWA=Libs/bwa/utils.o Libs/bwa/kthread.o Libs/bwa/kstring.o Libs/bwa/ksw.o Libs/bwa/bwt.o Libs/bwa/bntseq.o Libs/bwa/bwa.o Libs/bwa/bwamem.o Libs/bwa/bwamem_pair.o Libs/bwa/bwamem_extra.o Libs/bwa/malloc_wrap.o Libs/bwa/QSufSort.o Libs/bwa/bwt_gen.o Libs/bwa/rope.o Libs/bwa/rle.o Libs/bwa/is.o Libs/bwa/bwtindex.o Libs/bwa/bwashm.o Libs/bwa/bwase.o Libs/bwa/bwaseqio.o Libs/bwa/bwtgap.o Libs/bwa/bwtaln.o Libs/bwa/bamlite.o Libs/bwa/bwape.o Libs/bwa/kopen.o Libs/bwa/pemerge.o Libs/bwa/maxk.o Libs/bwa/bwtsw2_core.o Libs/bwa/bwtsw2_main.o Libs/bwa/bwtsw2_aux.o Libs/bwa/bwt_lite.o Libs/bwa/bwtsw2_chain.o Libs/bwa/fastmap.o Libs/bwa/bwtsw2_pair.o

TARGETS=srnaMapper # srnaCollapser srnaBuilder

all: $(TARGETS)

srnaCollapser: srnaCollapser.c
	$(CC) $(OFLAGS) -Wall -Wextra -o srnaCollapser srnaCollapser.c

srnaMapper: srnaMapper.c
	$(CC) $(OFLAGS) -Wall -Wextra -o srnaMapper.o -c srnaMapper.c
	$(CC) $(OFLAGS) -Wall -Wextra -o srnaMapper srnaMapper.o $(BWA) -lm -lz -lpthread -lrt


srnaBuilder: srnaBuilder.cpp
	$(CPPC) $(CFLAGS) -o srnaBuilder srnaBuilder.cpp $(LFLAGS)

clean:
	rm -f $(TARGETS)
