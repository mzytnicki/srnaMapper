#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include <locale.h>

#include "Libs/bwa/bwt.h"
#include "Libs/bwa/bwa.h"
#include "Libs/bwa/bwase.h"

#define DEPTH_SHORT_CUT 21

#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

#define N_NUCLEOTIDES    4
#define NUCLEOTIDES_BITS 2
#define NUCLEOTIDE_MASK  3
#define TRIPLET          3
#define N_TRIPLETS      64
#define TRIPLET_MASK    15

// The tree starts at the list of the TREE_BASE_SIZE-mers, which are N_TREE_BASE different k-mers.
#define TREE_BASE_SIZE   8
#define N_TREE_BASE  65536

#define N_STATES        0x100000
#define MANY_STATES     0x1000

#define PREPROCESSED_DEPTH 15

#define INIT_N_CELLS     0x1000000
#define INIT_N_QUALITIES 0x100000

#define NO_DATA 0
#define NO_QUALITY ((unsigned int) (-1))

#define BACKTRACE_SIZE    2
#define MATCH             0
#define MISMATCH          1
#define INSERTION         2
#define DELETION          3

#define CIGAR_SECONDARY_HIT 0x100
#define CIGAR_REVERSE 0x10

#define EDGE_STR_LENGTH 32
#define EDGE_SEQ_LENGTH 26
#define EDGE_LENGTH_LENGTH 6
#define MAX_EDGE_LENGTH 13
#define MAX_SW_LENGTH 20
#define MAX_SW_COST_SIZE 5
#define MAX_SW_N_STATES 10
#define MAX_SW_N_SEQUENCES 5
#define SW_BACKTRACE_SIZE_SIZE 5
#define SW_MATCH         0
#define SW_MISMATCH      1
#define SW_INSERTION     2
#define SW_DELETION      3

#define MAX_CIGAR_SIZE 510
#define MAX_QNAME_SIZE 255

// Should not do this, but whatever... Not really used...
#define MAX_READ_LENGTH 500

#define BWT_BUFFER_SIZE 10000

typedef unsigned long count_t;

const char DNA5_TO_CHAR [5] = { 'A', 'C', 'G', 'T', 'N' };
const unsigned char DNA5_TO_INT_REV [2][5] = { { 0, 1, 2, 3, 4 }, { 3, 2, 1, 0, 4 }};
const char *DNA5_TO_CHAR_REV[] = {"ACGTN", "TGCAN"};

const int  CHAR_TO_DNA5 [127] = {
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //   0
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //  10
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //  20
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //  30
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //  40
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //  50
 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, //  60
 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, //  70
 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, //  80
 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, //  90
 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, // 100
 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, // 110
 4, 4, 4, 4, 4, 4, 4           // 120
};
const char COMP [127] = {
 0,   0, 0,   0,   0,   0,   0,   0, 0,   0, //   0
 0,   0, 0,   0,   0,   0,   0,   0, 0,   0, //  10
 0,   0, 0,   0,   0,   0,   0,   0, 0,   0, //  20
 0,   0, 0,   0,   0,   0,   0,   0, 0,   0, //  30
 0,   0, 0,   0,   0,   0,   0,   0, 0,   0, //  40
 0,   0, 0,   0,   0,   0,   0,   0, 0,   0, //  50
 0,   0, 0,   0,   0, 'T',   0, 'G', 0,   0, //  60
 0, 'C', 0,   0,   0,   0,   0,   0, 0,   0, //  70
 0,   0, 0,   0, 'A', 'A',   0,   0, 0,   0, //  80
 0,   0, 0,   0,   0,   0,   0, 'T', 0, 'G', //  90
 0,   0, 0, 'G',   0,   0,   0,   0, 0,   0, // 100
 0,   0, 0,   0,   0, 'A', 'A',   0, 0,   0, // 110
 0,   0, 0,   0,   0,   0,   0               // 120
};

const char CIGAR [4] = { '=', 'X', 'I', 'D' };

bwt_t *bwt;
uint8_t *pac;
bntseq_t *bns;
unsigned long nReads;

typedef struct {
  bwtint_t k, l;
} bwtinterval_t;


#endif
