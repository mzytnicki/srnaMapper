// ==========================================================================
//                              srnaMapper
// ==========================================================================
// Copyright (C) 2019 Matthias Zytnicki, INRA
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Matthias Zytnicki or the INRA nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL MATTHIAS ZYTNICKI NOTR THE INRA BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Matthias Zytnicki <matthias.zytnicki@inra.fr>
// ==========================================================================

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

#define N_STATES        0x3000
#define MANY_STATES     0x1000

#define PREPROCESSED_DEPTH 15

#define INIT_N_CELLS     0x1000000
#define INIT_N_QUALITIES 0x100000

#define NO_DATA 0
#define NO_QUALITY ((unsigned int) (-1))

#define MATCH             0
#define MISMATCH          4
#define INSERTION         8
#define DELETION         12
#define BACKTRACE_OFFSET  2
#define BACKTRACE_MASK   12
#define PREPROCESSED     16

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
#define SW_BACKTRACE_SIZE 2
#define SW_BACKTRACE_SIZE_SIZE 5
#define SW_MATCH         0
#define SW_MISMATCH      1
#define SW_INSERTION     2
#define SW_DELETION      3

#define MAX_CIGAR_SIZE 510
#define MAX_QNAME_SIZE 255

#define BWT_BUFFER_SIZE 10

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


void trimSequence (size_t l, char *s) {
  s[l-1] = 0;
}

char *reverseSequence (char *sequence, char *revSequence, size_t size) {
  assert(strlen(sequence) == size);
  for (size_t i = 0; i < size; ++i) {
    revSequence[i] = sequence[size-i-1];
  }
  revSequence[size] = 0;
  return revSequence;
}

char *reverseComplementSequence (char *sequence, char *revSequence, size_t size) {
  for (size_t i = 0; i < size; ++i) {
    revSequence[i] = COMP[(int) sequence[size-i-1]];
  }
  revSequence[size] = 0;
  return revSequence;
}

unsigned short getNucleotide (unsigned long sequence, size_t position) {
  return ((sequence >> (position * NUCLEOTIDES_BITS)) & NUCLEOTIDE_MASK);
}

void printSequence (uint64_t sequence, size_t length) {
  for (unsigned int i = 0; i < length; ++i) {
    printf("%c", "ACGT"[sequence & NUCLEOTIDE_MASK]);
    sequence >>= NUCLEOTIDES_BITS;
  }
}

void printSequenceLong (unsigned short *sequence, size_t length) {
  for (unsigned int i = 0; i < length; ++i) {
    printf("%c", "ACGT"[sequence[i]]);
  }
}

typedef struct {
  unsigned long int nReads;
  unsigned long int nShortReads;
  unsigned long int nDown;
  unsigned long int nDownPreprocessed;
  unsigned long int nBufferCalls;
  unsigned long int nBufferCallSucesses;
  size_t            maxNStates;
} stats_t;

stats_t *stats;

void initializeStats () {
  stats->nReads              = 0;
  stats->nShortReads         = 0;
  stats->nDown               = 0;
  stats->nDownPreprocessed   = 0;
  stats->nBufferCalls        = 0;
  stats->nBufferCallSucesses = 0;
  stats->maxNStates          = 0;
}

void printStats () {
  char *savedLocale;
  savedLocale = setlocale (LC_ALL, NULL);
  setlocale(LC_NUMERIC, "");
  printf("Very small sequences: %'lu/%'lu\n", stats->nShortReads, stats->nReads);
  printf("# max states %'zu/%'i\n", stats->maxNStates, N_STATES);
  printf("# buffer call successes %'lu/%'lu (%i%%)\n", stats->nBufferCallSucesses, stats->nBufferCalls, (stats->nBufferCalls == 0)? 0: (int) (round(((double) stats->nBufferCallSucesses) / stats->nBufferCalls)));
  setlocale(LC_ALL, savedLocale);
}

typedef struct {
    char *readsFileNames[255];
    char *genomeFileName;
    char *outputReadsFileName;
    char *outputSamFileName;
    unsigned int nReadsFiles;
    size_t       maxNErrors;
    unsigned int lowComplexityThreshold;
    unsigned int maxNHits;
} parameters_t;

parameters_t *parameters;

void printUsage () {
  puts("srnaCollapser [-h] -r reads -g genome -o filename [-c filename] [-f filter] [-e #errors] [-n #max_hits]");
}

int parseCommandLine (int argc, char const **argv) {
  char *endptr;
  parameters->outputReadsFileName    = NULL;
  parameters->outputSamFileName      = NULL;
  parameters->nReadsFiles            = 0;
  parameters->maxNErrors             = 2;
  parameters->lowComplexityThreshold = 6;
  parameters->maxNHits               = 5;
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-h") == 0) {
      printUsage();
      return EXIT_FAILURE;
    }
    if (strcmp(argv[i], "-r") == 0) {
      ++i;
      parameters->readsFileNames[parameters->nReadsFiles] = strdup(argv[i]);
      ++parameters->nReadsFiles;
    }
    else if (strcmp(argv[i], "-g") == 0) {
      ++i;
      parameters->genomeFileName = strdup(argv[i]);
    }
    else if (strcmp(argv[i], "-o") == 0) {
      ++i;
      parameters->outputSamFileName = strdup(argv[i]);
    }
    else if (strcmp(argv[i], "-c") == 0) {
      ++i;
      parameters->outputReadsFileName = strdup(argv[i]);
    }
    else if (strcmp(argv[i], "-f") == 0) {
      ++i;
      parameters->lowComplexityThreshold = strtol(argv[i], &endptr, 10);
    }
    else if (strcmp(argv[i], "-e") == 0) {
      ++i;
      parameters->maxNErrors = strtol(argv[i], &endptr, 10);
    }
    else if (strcmp(argv[i], "-n") == 0) {
      ++i;
      parameters->maxNHits = strtol(argv[i], &endptr, 10);
    }
    else {
      printf("Cannot understand parameter '%s'\n", argv[i]);
      printUsage();
      return EXIT_FAILURE;
    }
  }
  if (parameters->nReadsFiles == 0) {
    printUsage();
    puts("Reads file is missing.\nExiting.");
    return EXIT_FAILURE;
  }
  if (parameters->genomeFileName == NULL) {
    printUsage();
    puts("Genome index file is missing.\nExiting.");
    return EXIT_FAILURE;
  }
  if (parameters->outputSamFileName == NULL) {
    printUsage();
    puts("Output SAM file is missing.\nExiting.");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

typedef struct {
  FILE *file;
  char *qname;
  unsigned int backtraceSize;
  char *backtraceCigar;
  unsigned char *backtraceLengths;
  char *forwardCigar;
  char *backwardCigar;
  char *forwardSeq;
  char *backwardSeq;
  char *forwardQual;
  char *backwardQual;
  bool isBackwardSet;
} outputSam_t;

void createOutputSam (outputSam_t *outputSam, unsigned int size) {
  outputSam->isBackwardSet    = false;
  outputSam->qname            = (char *) malloc(MAX_QNAME_SIZE * sizeof(char));
  outputSam->backtraceCigar   = (char *) malloc(MAX_CIGAR_SIZE * sizeof(char));
  outputSam->backtraceLengths = (unsigned char *) malloc(MAX_CIGAR_SIZE * sizeof(unsigned char));
  outputSam->forwardCigar     = (char *) malloc(MAX_CIGAR_SIZE * sizeof(char));
  outputSam->backwardCigar    = (char *) malloc(MAX_CIGAR_SIZE * sizeof(char));
  outputSam->forwardSeq       = (char *) malloc((size+1) * sizeof(char));
  outputSam->backwardSeq      = (char *) malloc((size+1) * sizeof(char));
  outputSam->forwardQual      = (char *) malloc((size+1) * sizeof(char));
  outputSam->backwardQual     = (char *) malloc((size+1) * sizeof(char));
}

void computeCigar (outputSam_t *outputSam) {
  size_t forwardCigarSize = 0;
  size_t backwardCigarSize = 0;
  //TODO Optimize this
  forwardCigarSize = backwardCigarSize = 0;
  for (size_t backtraceId = 0; backtraceId < outputSam->backtraceSize; ++backtraceId) {
    forwardCigarSize  += sprintf(outputSam->forwardCigar+forwardCigarSize,   "%i%c", outputSam->backtraceLengths[backtraceId], outputSam->backtraceCigar[backtraceId]);
    backwardCigarSize += sprintf(outputSam->backwardCigar+backwardCigarSize, "%i%c", outputSam->backtraceLengths[outputSam->backtraceSize-backtraceId-1], outputSam->backtraceCigar[outputSam->backtraceSize-backtraceId-1]);
  }
  outputSam->forwardCigar[forwardCigarSize] = 0;
  outputSam->backwardCigar[backwardCigarSize] = 0;
}

void setCigarNoError (outputSam_t *outputSam, unsigned int readLength) {
  sprintf(outputSam->forwardCigar, "%uM", readLength);
  sprintf(outputSam->backwardCigar, "%uM", readLength);
}

void computeReverseComplement (outputSam_t *outputSam) {
  size_t readLength = strlen(outputSam->forwardSeq);
  outputSam->isBackwardSet = true;
  reverseComplementSequence(outputSam->forwardSeq, outputSam->backwardSeq, readLength);
  reverseSequence(outputSam->forwardQual, outputSam->backwardQual, readLength);
}

void freeOutputSam (outputSam_t *outputSam) {
  free(outputSam->qname);
  free(outputSam->forwardCigar);
  free(outputSam->backwardCigar);
  free(outputSam->forwardSeq);
  free(outputSam->backwardSeq);
  free(outputSam->forwardQual);
  free(outputSam->backwardQual);
}

/******* Edge type *******/
/**
 * An edge is the link from one cell to one of its children.
 * It also contains the sequence for a cell to the child.
 * The sequence is stored in a int, with the first nucleotides being the sequence, and the last ones the size
 * It is:
 *   - the sequence: mix of sequence and length
 *   - the id of the next cell
 */
typedef uint32_t sequence_t;

typedef struct {
  sequence_t sequence: EDGE_SEQ_LENGTH;
  sequence_t length:   EDGE_LENGTH_LENGTH;
  uint64_t cellId;
} edge_t;

void createEdge (edge_t *edge) {
  edge->length   = 0;
  edge->sequence = 0;
  edge->cellId = NO_DATA;
}

bool isSetEdge (edge_t *edge) {
  return (edge->length > 0);
}

bool isFullEdge (edge_t *edge) {
  return (edge->length == MAX_EDGE_LENGTH);
}

void setEdge (edge_t *edge, sequence_t sequence, sequence_t length, uint64_t cellId) {
  edge->sequence = sequence;
  edge->length   = length;
  edge->cellId   = cellId;
}

void unsetEdge(edge_t *edge) {
  edge->sequence = 0;
  edge->length   = 0;
  edge->cellId   = NO_DATA;
}

void setCellIdEdge (edge_t *edge, uint64_t cellId) {
  edge->cellId = cellId;
}

void addEdgeNucleotide (edge_t *edge, unsigned short nucleotide) {
  assert(! isFullEdge(edge));
  //printf("    adding nucl. %u to %i @ %i\n", nucleotide, edge->sequence, edge->length);
  edge->sequence |= nucleotide << (NUCLEOTIDES_BITS * edge->length);
  //printf("      now: %i\n", edge->sequence);
  ++edge->length;
}

unsigned short getEdgeNucleotide (edge_t *edge, sequence_t length) {
  assert(length <= edge->length);
  assert(length <= MAX_EDGE_LENGTH);
  return getNucleotide(edge->sequence, length);
}

void edgeRemoveFirst (edge_t *edge, size_t length) {
  assert(length <= edge->length);
  edge->length -= length;
  edge->sequence >>= (length * NUCLEOTIDES_BITS);
}

void edgeRemoveLast (edge_t *edge, size_t newSize) {
  assert(newSize <= edge->length);
  edge->length = newSize;
  edge->sequence &= ((1 << (newSize * NUCLEOTIDES_BITS)) - 1);
}

void edgeSetCellId (edge_t *edge, uint64_t cellId) {
  edge->cellId = cellId;
}

void printEdge (edge_t *edge) {
  if (! isSetEdge(edge)) {
    printf("-X->");
    return;
  }
  printf("-");
  printSequence(edge->sequence, edge->length);
  printf("-> %" PRIu64, edge->cellId);
}

/******* Cell type *******/
/**
 * A cell is a prefix of the reads tree.
 * It is:
 *   - 4 children: an array of 4 indices of the children array in the tree structure
 *   - the read counts: an array of int
 */

typedef struct {
  edge_t edges [N_NUCLEOTIDES];
  count_t *counts;
} cell_t;

void createCell (cell_t *cell) {
  for (unsigned short i = 0; i < N_NUCLEOTIDES; ++i) {
    createEdge(&cell->edges[i]);
  }
  //printf("%u %zu %p\n", parameters->nReadsFiles, sizeof(count_t), cell);
  cell->counts = (count_t *) calloc(parameters->nReadsFiles, sizeof(count_t));
}

void freeCell (cell_t *cell) {
  free(cell->counts);
}

void addEdge (cell_t *cell, sequence_t sequence, sequence_t length, uint64_t childId) {
  unsigned short nucleotide = sequence & NUCLEOTIDE_MASK;
  setEdge(&cell->edges[nucleotide], sequence, length, childId);
}

edge_t *getFirstEdge (cell_t *cell) {
  edge_t *edge;
  for (unsigned short edgeId = 0; edgeId < N_NUCLEOTIDES; ++edgeId) {
    edge = &cell->edges[edgeId];
    if (isSetEdge(edge)) {
      return edge;
    }
  }
  return NULL;
}

/**
 * Split an edge into two.  Add the remaining part of the edge to a new cell.
 * When newEdgeId == N_NUCLEOTIDES, the nucleotide at the split point is not known
 */
void splitEdge (edge_t *edge, cell_t *newCell, uint64_t newCellId, sequence_t length, unsigned short newEdgeId) {
  edge_t *newEdge;
  //printf("   Splitting edge %p %i into %u -> %" PRIu64 " at point %zu and id %i\n", edge, edge->sequence, edge->length, edge->cellId, length, newEdgeId);
  assert(edge->cellId != NO_DATA);
  assert(length < edge->length);
  assert(newCellId != NO_DATA);
  if (newEdgeId == N_NUCLEOTIDES) newEdgeId = getEdgeNucleotide(edge, length);
  newEdge = &newCell->edges[newEdgeId];
  newEdge->sequence = edge->sequence;
  newEdge->length = edge->length;
  edgeRemoveFirst(newEdge, length);
  edgeRemoveLast(edge, length);
  newEdge->cellId = edge->cellId;
  edge->cellId = newCellId;
  //printf("   split edge %p into %i (%u) -> %" PRIu64 " and %p %i (%u) -> %" PRIu64 "\n", edge, edge->sequence, edge->length, edge->cellId, newEdge, newEdge->sequence, newEdge->length, newEdge->cellId);
}

unsigned short getNChildren (cell_t *cell) {
  unsigned short nChildren = 0;
  for (unsigned int edgeId = 0; edgeId < N_NUCLEOTIDES; ++edgeId) {
    if (isSetEdge(&cell->edges[edgeId])) {
      ++nChildren;
    }
  }
  return nChildren;
}

bool isTerminal (cell_t *cell) {
  for (unsigned int edgeId = 0; edgeId < N_NUCLEOTIDES; ++edgeId) {
    if (isSetEdge(&cell->edges[edgeId])) {
      return false;
    }
  }
  return true;
}

void printCell (cell_t *cell) {
  for (unsigned short i = 0; i < N_NUCLEOTIDES; ++i) {
    if (isSetEdge(&cell->edges[i])) {
      printEdge(&cell->edges[i]);
      printf("  ");
    }
  }
}


/******* Quality type *******/
/**
 * A quality the quality of a set of reads with the same sequence.
 * It is:
 *   - the qualities: an array of char*
 *   - the size of the array
 * The quality array is roughly of the same size as the tree array.
 * If the quality is set for a cell, the corresponding index should not be null.
 */

typedef struct {
  char       **qualities;
  //uint64_t    *cellIds;
  //unsigned int nQualities;
  unsigned int nAllocated;
} quality_t;

void createQualities (quality_t *qualities) {
  qualities->nAllocated = INIT_N_QUALITIES;
  qualities->qualities  = (char **) calloc(qualities->nAllocated, sizeof(char *));
  //qualities->nQualities = 0;
}

void freeQualities (quality_t *qualities) {
  for (unsigned int i = 0; i < qualities->nAllocated; ++i) {
    if (qualities->qualities[i] != NULL) {
      free(qualities->qualities[i]);
    }
  }
  free(qualities->qualities);
}

unsigned int findQualityId (const quality_t *qualities, uint64_t cellId) {
  if ((cellId >= qualities->nAllocated) || (qualities->qualities[cellId] == NULL)) return NO_QUALITY;
  return cellId;
}

char *findQuality (const quality_t *qualities, uint64_t cellId) {
  unsigned int qualityId = findQualityId(qualities, cellId);
  if (qualityId == NO_QUALITY) return NULL;
  //printf("\t\t\tFinding quality for %" PRIu64 ": %u\n", cellId, qualityId);
  return qualities->qualities[qualityId]; 
}

void addQuality (quality_t *qualities, uint64_t cellId, size_t l, char *quality) {
  //printf("Adding quality %s at %" PRIu64 ".\n", quality, cellId);
  //TODO: Is the "if" useful?
  size_t previousNAllocated;
  if (cellId >= qualities->nAllocated) {
    previousNAllocated = qualities->nAllocated;
    while (cellId >= qualities->nAllocated) {
      qualities->nAllocated *= 2;
    }
    if ((qualities->qualities = (char **) realloc(qualities->qualities, qualities->nAllocated * sizeof(char *))) == NULL) {
      printf("Cannot allocate memory for qualities of size %u.\nExiting.\n", qualities->nAllocated);
      exit(EXIT_FAILURE);
    }
    for (size_t i = previousNAllocated + 1; i < qualities->nAllocated; ++i) {
      qualities->qualities[i] = NULL;
    }
  }
  qualities->qualities[cellId] = strndup(quality, l);
}

void replaceQuality (quality_t *qualities, unsigned int qualityId, size_t l, char *quality) {
  //printf("Quality: %zu vs %zu\n", strlen(qualities->qualities[qualityId]), strlen(quality));
  assert(strlen(qualities->qualities[qualityId]) == strlen(quality));
  assert(strlen(quality) == l);
  for (size_t i = 0; i < l; ++i) {
    qualities->qualities[qualityId][i] = MAX(qualities->qualities[qualityId][i], quality[i]);
  }
}

void _setQuality (quality_t *qualities, uint64_t cellId, size_t l, char *quality) {
  unsigned qualityId = findQualityId(qualities, cellId);
  if (qualityId == NO_QUALITY) {
    addQuality(qualities, cellId, l, quality);
  }
  else {
    replaceQuality(qualities, qualityId, l, quality);
  }
}

/******* Tree type *******/
/**
 * A tree is stores all the prefixes of the reads.
 * It is:
 *   - a depth
 *   - an array of cells
 *   - the number of cells
 *   - the size of the array
 *   - an instance of the quality structure
 * The first N_TREE_BASE cells correspond to the all the combinations of TREE_BASE_SIZE-mers.
 * In other words, the 1-mers, 2-mers, ..., TREE_BASE_SIZE-1-mers are not represented.
 */

typedef struct {
  size_t    depth;
  cell_t   *cells;
  uint64_t  nCells;
  uint64_t  nAllocated;
  quality_t qualities;
} tree_t;

void createTree (tree_t *tree) {
  tree->nAllocated = INIT_N_CELLS;
  tree->cells = (cell_t *) malloc(tree->nAllocated * sizeof(cell_t));
  tree->depth = TREE_BASE_SIZE;
  tree->nCells = N_TREE_BASE;
  for (size_t i = 0; i < N_TREE_BASE; ++i) {
    createCell(&tree->cells[i]);
  }
  createQualities(&tree->qualities);
}

void _computeTreeStats (const tree_t *tree, unsigned int **stats, unsigned int *statsSum, unsigned int *branchSizes, unsigned int branchSize, cell_t *cell, size_t depth, unsigned int *nNodes, unsigned int *nQualities) {
  edge_t *edge;
  size_t length;
  unsigned int c = 0;
  for (short edgeId = 0; edgeId < N_NUCLEOTIDES; ++edgeId) {
    edge = &cell->edges[edgeId];
    if (edge->cellId != NO_DATA) {
      length = edge->length;
      ++c;
      ++(*nNodes);
      if (findQualityId(&tree->qualities, edge->cellId) != NO_QUALITY) {
        ++(*nQualities);
      }
    }
  }
  if (c == 1) {
    branchSize += length;
  }
  else {
    ++branchSizes[branchSize];
    branchSize = 0;
  }
  for (short edgeId = 0; edgeId < N_NUCLEOTIDES; ++edgeId) {
    edge = &cell->edges[edgeId];
    if (edge->cellId != NO_DATA) {
      _computeTreeStats(tree, stats, statsSum, branchSizes, branchSize, &tree->cells[edge->cellId], depth+edge->length, nNodes, nQualities);
    }
  }
  //printf("At depth %zu/%zu with %i nucleotides\n", depth, tree->depth, c); fflush(stdout);
  ++stats[depth][c];
  ++statsSum[c];
}

void computeTreeStats (const tree_t *tree) {
  //unsigned int **stats = (unsigned int **) malloc((tree->depth+1) * N_NUCLEOTIDES * sizeof(unsigned int *));
  unsigned int **stats = (unsigned int **) malloc((tree->depth+1) * sizeof(unsigned int *));
  unsigned int statsSum[N_NUCLEOTIDES] = { 0, 0, 0, 0 };
  unsigned int *branchSizes = (unsigned int *) calloc(tree->depth+1, sizeof(unsigned int));
  unsigned int nNodes = 0;
  unsigned int nQualities = 0;
  unsigned int s;
  for (size_t depth = 0; depth <= tree->depth; ++depth) {
    stats[depth] = (unsigned int *) calloc(N_NUCLEOTIDES + 1, sizeof(unsigned int));
  }
  for (size_t cellId = 0; cellId < N_TREE_BASE; ++cellId) {
    _computeTreeStats(tree, stats, statsSum, branchSizes, 0, &tree->cells[cellId], TREE_BASE_SIZE, &nNodes, &nQualities);
  }
  puts("Stats on tree");
  for (size_t depth = TREE_BASE_SIZE; depth <= tree->depth; ++depth) {
    printf("%zu:", depth);
    s = 0;
    for (size_t nChildren = 0; nChildren < N_NUCLEOTIDES; ++nChildren) {
      s += stats[depth][nChildren];
    }
    if (s != 0) {
      for (size_t nChildren = 0; nChildren < N_NUCLEOTIDES; ++nChildren) {
        printf("\t%zu:%u (%.f%%)", nChildren, stats[depth][nChildren], ((float) stats[depth][nChildren]) / s * 100);
      }
    }
    printf("\n");
  }
  printf("sum:");
  s = 0;
  for (size_t nChildren = 0; nChildren < N_NUCLEOTIDES; ++nChildren) {
    s += statsSum[nChildren];
  }
  for (size_t nChildren = 0; nChildren < N_NUCLEOTIDES; ++nChildren) {
    printf("\t%zu:%u (%.f%%)", nChildren, statsSum[nChildren], ((float) statsSum[nChildren]) / s * 100);
  }
  printf("\n");
  puts("Stats on branch sizes");
  for (size_t depth = 0; depth <= tree->depth; ++depth) {
    printf("%zu: %u\n", depth, branchSizes[depth]);
  }
  printf("%u nodes with %u qualities (%.f%%)\n", nNodes, nQualities, ((float) nQualities) / nNodes * 100);
  for (size_t depth = 0; depth <= tree->depth; ++depth) {
    free(stats[depth]);
  }
  free(stats);
  free(branchSizes);
}

void freeTree (tree_t *tree) {
  tree->nAllocated = 0;
  tree->depth = 0;
  tree->nCells = 0;
  for (uint64_t i = 0; i < tree->nCells; ++i) {
    freeCell(&tree->cells[i]);
  }
  free(tree->cells);
  freeQualities(&tree->qualities);
}

void setQuality (tree_t *tree, size_t cellId, size_t l, char *quality, unsigned int fileId) {
  ++tree->cells[cellId].counts[fileId];
  _setQuality(&tree->qualities, cellId, l, quality);
}

uint64_t addCell (tree_t *tree) {
  if (tree->nCells == tree->nAllocated) {
    tree->nAllocated *= 2;
    //printf("reallocating cells...\n");
    if ((tree->cells = (cell_t *) realloc(tree->cells, tree->nAllocated * sizeof(cell_t))) == NULL) {
      printf("Cannot allocate memory for tree of size %lu.\nExiting.\n", tree->nAllocated);
      exit(EXIT_FAILURE);
    }
  }
  //printf("ncells: %zu / %zu\n", tree->nCells, tree->nAllocated);
  createCell(&tree->cells[tree->nCells]);
  ++tree->nCells;
  return tree->nCells-1;
}

/**
 * Create a new node, split an edge into two, and add the remaining part of the edge to then new cell.
 */
uint64_t splitEdgeTree (tree_t *tree, edge_t *edge, sequence_t length, unsigned short newEdgeId) {
  uint64_t newCellId = addCell(tree);
  splitEdge(edge, &tree->cells[newCellId], newCellId, length, newEdgeId);
  return newCellId;
}

bool isCellUnbranched (const tree_t *tree, cell_t *cell) {
  bool foundOneChild = false;
  unsigned short child = 0;
  //printf("\t\tCell unbranched: ");
  //printCell(cell);
  //printf("\n");
  for (unsigned int edgeId = 0; edgeId < N_NUCLEOTIDES; ++edgeId) {
    if (isSetEdge(&cell->edges[edgeId])) {
      //printEdge(&cell->edges[edgeId]);
      if (foundOneChild) {
        //printf("\t\t\tNo\n");
        return false;
      }
      foundOneChild = true;
      child = edgeId;
    }
  }
  if (foundOneChild) {
    return isCellUnbranched(tree, &tree->cells[cell->edges[child].cellId]);
  }
  //printf("\t\t\tYes\n");
  return true;
}

/**
 * Add the rest of the sequence and extend the edge
 */
uint64_t addSequenceAdd (tree_t *tree, uint64_t cellId, char *sequence, int sequenceId, edge_t *edge) {
  unsigned short nucleotide;
  //printf("Adding sequence %s @%" PRIu64 "\n", sequence, cellId);
  for (; sequenceId >= 0; --sequenceId) {
    nucleotide = CHAR_TO_DNA5[(int) sequence[sequenceId]];
    if (edge == NULL) {
      edge = &tree->cells[cellId].edges[nucleotide];
    }
    //printf("  seq id: %d, edge: %p (%i), nucleotide: %c\n", sequenceId, edge, edge->length, DNA5_TO_CHAR[nucleotide]);
    addEdgeNucleotide(edge, nucleotide);
    if (edge->length == MAX_EDGE_LENGTH) {
      cellId = addCell(tree);
      //printf("  to the end: %" PRIu64 "\n", cellId);
      setCellIdEdge(edge, cellId);
      edge = NULL;
    }
  }
  // set node
  if (edge != NULL) {
      cellId = addCell(tree);
      setCellIdEdge(edge, cellId);
  }
  //printf("  over with %" PRIu64 "\n", cellId);
  return cellId;
  /*
  cell_t *cell = &tree->cells[cellId];
  uint64_t newCellId = cell->children[childId];
  if (newCellId != NO_DATA) {
    return newCellId;
  }
  newCellId = addCell(tree);
  cell = &tree->cells[cellId]; // may be reallocated!
  cell->children[childId] = newCellId;
  return newCellId;
  */
}


uint64_t addSequenceFollow (tree_t *tree, uint64_t cellId, char *sequence, int sequenceId) {
  size_t edgeLength = 0;
  unsigned short sequenceNucleotide, edgeNucleotide;
  edge_t *edge = NULL;
  //printf("Following sequence %s @%" PRIu64 "\n", sequence, cellId); fflush(stdout);
  for (; sequenceId >= 0; --sequenceId) {
    sequenceNucleotide = CHAR_TO_DNA5[(int) sequence[sequenceId]];
    if (edge == NULL) {
      //printf("  new edge\n"); fflush(stdout);
      edge = &tree->cells[cellId].edges[sequenceNucleotide];
    }
    //printf("  seq id: %d, edge: %p (%i/%zu) -> %" PRIu64 ", nucleotide: %c\n", sequenceId, edge, edge->length, edgeLength, edge->cellId, DNA5_TO_CHAR[sequenceNucleotide]); fflush(stdout);
    if (edge->cellId == NO_DATA) {
      return addSequenceAdd(tree, cellId, sequence, sequenceId, edge);
    }
    edgeNucleotide = getEdgeNucleotide(edge, edgeLength); fflush(stdout);
    //printf("  edge %p, value: %i, len: %zu/%u, nucleotide: %c (%u)\n", edge, edge->sequence, edgeLength, edge->length, DNA5_TO_CHAR[edgeNucleotide], edgeNucleotide); fflush(stdout);
    if (sequenceNucleotide == edgeNucleotide) {
      //printf("  following to %" PRIu64 "\n", edge->cellId); fflush(stdout);
      ++edgeLength;
      if (edgeLength == edge->length) {
        cellId = edge->cellId;
        //printf("  to the end -> %" PRIu64 "\n", cellId); fflush(stdout);
        edge = NULL;
        edgeLength = 0;
      }
    }
    else {
      //printf("  split @ %zu\n", edgeLength); fflush(stdout);
      cellId = edge->cellId;
      if (edgeLength == 0) {
        edge = &tree->cells[cellId].edges[sequenceNucleotide];
        //printf("  edge is %p\n", edge); fflush(stdout);
        return addSequenceAdd(tree, cellId, sequence, sequenceId, edge);
      }
      else {
        cellId = splitEdgeTree(tree, edge, edgeLength, edgeNucleotide);
        edge   = &tree->cells[cellId].edges[sequenceNucleotide];
        return addSequenceAdd(tree, cellId, sequence, sequenceId, edge);
      }
    }
  }
  // set node
  //printf("  over with %zu\n", edgeLength); fflush(stdout);
  if (edgeLength != 0) {
    cellId = splitEdgeTree(tree, edge, edgeLength, N_NUCLEOTIDES);
  }
  return cellId;
}

bool addSequence (tree_t *tree, size_t l, char *sequence, char *quality, unsigned int fileId) {
  uint64_t cellId = 0;
  int sequenceId = l - 1;
  assert(strlen(sequence) == strlen(quality));
  assert(strlen(quality) == l);
  if (l < TREE_BASE_SIZE) {
    return false;
  }
  for (int i = 0; i < TREE_BASE_SIZE; ++i, --sequenceId) {
    cellId <<= NUCLEOTIDES_BITS;
    cellId += CHAR_TO_DNA5[(int) sequence[sequenceId]];
  }
  //printf("First id: %lu, %s\n", cellId, sequence);
  assert(cellId < N_TREE_BASE);
  cellId = addSequenceFollow(tree, cellId, sequence, sequenceId);
  setQuality(tree, cellId, l, quality, fileId);
  tree->depth = MAX(tree->depth, l);
  return true;
}

/**
 * Print the last nucleotides (close to the leaves) of the tree.
 */
void __printTree (const tree_t *tree, FILE *outFile, uint64_t *readId, char *read, size_t readPos, uint64_t cellId) {
  uint64_t nextCellId;
  edge_t *edge;
  cell_t *cell = &tree->cells[cellId];
  char *quality;
  //printf("Read: %s at %"PRIu64 "\n", read+tree->depth-readPos, cellId);
  if ((quality = findQuality(&tree->qualities, cellId)) != NULL) {
    //printf("\tGot quality\n"); fflush(stdout);
    ++(*readId);
    fprintf(outFile, "@read%" PRIu64 "_", *readId);
    for (unsigned int fileId = 0; fileId < parameters->nReadsFiles; ++fileId) {
      fprintf(outFile, "_%lu", cell->counts[fileId]);
    }
    fprintf(outFile, "\n%s\n+\n%s\n", read+tree->depth-readPos, quality);
    //assert(strlen(read+tree->depth-readPos) == strlen(cell->quality));
  }
  for (unsigned short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    edge = &cell->edges[nucleotide];
    nextCellId = edge->cellId;
    if (nextCellId != NO_DATA) {
      for (size_t i = 0; i < edge->length; ++i) {
        assert(1 + readPos + i <= tree->depth);
        read[tree->depth-readPos-1-i] = DNA5_TO_CHAR[getEdgeNucleotide(edge, i)];
      }
      __printTree(tree, outFile, readId, read, readPos + edge->length, nextCellId);
    }
  }
}

/**
 * Print the first nucleotides (close to the root) of the tree, then call __printTree.
 */
void _printTree (const tree_t *tree, FILE *outFile, uint64_t *readId, char *read, size_t readPos, uint64_t cellId) {
  if (readPos == TREE_BASE_SIZE) {
    __printTree (tree, outFile, readId, read, readPos, cellId);
    return;
  }
  assert(cellId < N_TREE_BASE);
  cellId <<= NUCLEOTIDES_BITS;
  for (unsigned short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    read[tree->depth-readPos-1] = DNA5_TO_CHAR[nucleotide];
    _printTree(tree, outFile, readId, read, readPos+1, cellId+nucleotide);
  }
}

/**
 * Open/close file, allocate the memory, and call _printTree
 */
int printTree (char *fileName, const tree_t *tree) {
  FILE *outFile = fopen(fileName, "w");
  if (outFile == NULL) return EXIT_FAILURE;
  uint64_t readId = 0;
  char *read = (char *) malloc((tree->depth+1) * sizeof(char));
  read[tree->depth] = 0;
  _printTree(tree, outFile, &readId, read, 0, 0);
  free(read);
  fclose(outFile);
  printf("Done with print tree.\n");
  return EXIT_SUCCESS;
}

int readReadsFile (char *fileName, tree_t *tree, unsigned int fileId) {
  FILE *inFile;
  char *line = NULL;
  char *sequence = NULL;
  char *quality = NULL;
  size_t len = 0;
  ssize_t nRead;
  inFile = fopen(fileName, "r");
  if (inFile == NULL) return EXIT_FAILURE;
  while ((nRead = getline(&line, &len, inFile)) != -1) {
    nRead = getline(&sequence, &len, inFile);
    if (nRead == -1) {
      fprintf(stderr, "Input file '%s' is corrupted.\nAborting.\n", fileName);
      return EXIT_FAILURE;
    }
    nRead = getline(&line, &len, inFile);
    if (nRead == -1) {
      fprintf(stderr, "Input file '%s' is corrupted.\nAborting.\n", fileName);
      return EXIT_FAILURE;
    }
    nRead = getline(&quality, &len, inFile);
    if (nRead == -1) {
      fprintf(stderr, "Input file '%s' is corrupted.\nAborting.\n", fileName);
      return EXIT_FAILURE;
    }
    assert(strlen(sequence) == strlen(quality));
    assert(strlen(sequence) == (unsigned long) nRead);
    trimSequence(nRead, sequence);
    trimSequence(nRead, quality);
    ++stats->nReads;
    if (! addSequence(tree, nRead-1, sequence, quality, fileId)) {
      ++stats->nShortReads;
    }
  }
  free(line);
  free(sequence);
  free(quality);
  fclose(inFile);
  return EXIT_SUCCESS;
}

/**
 * Count the 3-mer, and possibly filter
 */
bool __filterTree (const tree_t *tree, size_t readPos, uint64_t cellId, unsigned short prevTriplet, count_t *prevTripletCount) {
  edge_t *edge;
  uint64_t nextCellId;
  bool foundRead = false;
  bool thresholdReached = false;
  cell_t *cell = &tree->cells[cellId];
  count_t tripletCount [N_TRIPLETS];
  unsigned short triplet;
  size_t prevReadPos = readPos;
  unsigned short nucleotide;
  //printf("\tCurrent triplet: %u @ readPos: %zu\n", prevTriplet, readPos);
  if (findQuality(&tree->qualities, cellId) != NULL) {
    foundRead = true;
  }
  for (unsigned short edgeId = 0; edgeId < N_NUCLEOTIDES; ++edgeId) {
    edge = &cell->edges[edgeId];
    nextCellId = edge->cellId;
    //printf("\t\tTrying nucleotide %i\n", edgeId);
    if (nextCellId != NO_DATA) {
      thresholdReached = false;
      triplet = prevTriplet;
      memcpy(tripletCount, prevTripletCount, N_TRIPLETS * sizeof(count_t));
      readPos = prevReadPos;
      for (size_t edgeLength = 0; (edgeLength < edge->length) && (! thresholdReached); ++edgeLength, ++readPos) {
        nucleotide = getEdgeNucleotide(edge, edgeLength);
        triplet = ((triplet & TRIPLET_MASK) << NUCLEOTIDES_BITS) | nucleotide;
        //printf("\t\t\ttriplet: %u @ edgelen: %zu, nucleotide: %i\n", triplet, edgeLength, nucleotide);
        if (readPos >= TRIPLET-1) tripletCount[triplet] += 1;
        if (tripletCount[triplet] > parameters->lowComplexityThreshold) {
          //printf("\t\t\t\tDirect filtering\n");
          thresholdReached = true;
        }
     }
     if (thresholdReached) {
        unsetEdge(edge);
      }
      else {
        if (! __filterTree(tree, readPos+1, nextCellId, triplet, tripletCount)) {
          //printf("\t\tSecond filtering @ %zu\n", readPos);
          unsetEdge(edge);
        }
        else {
          foundRead = true;
        }
      }
    }
  }
  return foundRead;
}

/**
 * Initialize counts and run _filterTree
 */
bool _filterTree (const tree_t *tree, size_t readPos, uint64_t cellId, unsigned short triplet, count_t *tripletCount) {
  if (readPos == TREE_BASE_SIZE) {
    return __filterTree(tree, readPos, cellId, triplet, tripletCount);
  }
  unsigned short nextTriplet;
  bool foundRead = false;
  cellId <<= NUCLEOTIDES_BITS;
  triplet &= TRIPLET_MASK;
  triplet <<= NUCLEOTIDES_BITS;
  for (unsigned short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    nextTriplet = triplet | nucleotide;
    //printf("\tFirst filter @ %zu with nt %i and triplet %i\n", readPos, nucleotide, triplet);
    if (readPos >= TRIPLET-1) tripletCount[nextTriplet] += 1;
    if (_filterTree(tree, readPos+1, cellId | nucleotide, nextTriplet, tripletCount)) {
      foundRead = true;
    }
    if (readPos >= TRIPLET-1) tripletCount[nextTriplet] -= 1;
  }
  return foundRead;
}

/**
 * Initialize counts and run _filterTree
 */
unsigned int filterTree (const tree_t *tree) {
  count_t tripletCount [N_TRIPLETS];
  for (unsigned int i = 0; i < N_TRIPLETS; ++i) tripletCount[i] = 0;
  return _filterTree(tree, 0, 0, 0, tripletCount);
}

bwaidx_t *loadGenomeFile (char *indexName) {
  bwaidx_t *idx = bwa_idx_load(indexName, BWA_IDX_ALL);
  if (idx == NULL) {
    fprintf(stderr, "Index load failed.\n");
  }
  return idx;
}

/******* BWT-interval type *******/
/**
 * A BWT-interval coresponds to a suffix of the genome.
 * It is:
 *   - the index of the first interval of the SA
 *   - the index of the last+1 interval of the SA
 *   - an array of cells
 *   - the number of cells
 *   - the size of the array
 *   - an instance of the quality structure
 * The first N_TREE_BASE cells correspond to the all the combinations of TREE_BASE_SIZE-mers.
 * In other words, the 1-mers, 2-mers, ..., TREE_BASE_SIZE-1-mers are not represented.
 */

typedef struct {
  bwtint_t k, l;
} bwtinterval_t;

bool compareBwtIntervals (bwtinterval_t *i1, bwtinterval_t *i2) {
  return ((i1->k == i2->k) && (i1->l == i2->l));
}

//static unsigned long nIntervals;

/******* State type *******/
/**
 * A state is an prefix of the genome tree, with additional information to perform backtrace.
 * It is:
 *   - a BWT-interval
 *   - the current operation of the backtrace
 *   - the parent state for a backtrace
 */

typedef struct state_t state_t;

struct state_t {
  bwtinterval_t interval;
  unsigned char trace;
  struct state_t *previousState;
};

void setState (state_t *state, bwtinterval_t *interval, unsigned char trace, state_t *previousState) {
  state->interval      = *interval;
  state->trace         = trace;
  state->previousState = previousState;
}

bwtinterval_t *getStateInterval (state_t *state) {
  return &state->interval;
}

bool hasTrace (state_t *state, unsigned char traceType) {
  return ((state->trace & BACKTRACE_MASK) == traceType);
}

void printState (state_t *state, size_t maxDepth) {
  //assert(getStateInterval(state)->k <= getStateInterval(state)->l);
  size_t s = 0;
  char *tmpSeq1, *tmpSeq2;
  tmpSeq1 = (char *) malloc(255);
  tmpSeq2 = (char *) malloc(255);
  tmpSeq1[maxDepth] = 0;
  tmpSeq2[maxDepth] = 0;
  int is_rev;
  bwtint_t p1 = bwt_sa(bwt, getStateInterval(state)->k); // position on the forward-reverse coordinate
  p1 = bns_depos(bns, p1, &is_rev); // position on the forward strand; this may be the first base or the last base
  for (size_t i = 0; i < maxDepth; ++i) {
    bwtint_t p = is_rev? p1-i: p1+i;
    if (p == ULONG_MAX) {
      tmpSeq1[i] = 0;
      break;
    }
    tmpSeq1[i] = ((is_rev)? "TGCAN": "ACGTN")[_get_pac(pac, p)];
  }
  bwtint_t p2 = bwt_sa(bwt, getStateInterval(state)->l); // position on the forward-reverse coordinate
  p2 = bns_depos(bns, p2, &is_rev); // position on the forward strand; this may be the first base or the last base
  for (size_t i = 0; i < maxDepth; ++i) {
    bwtint_t p = is_rev? p2-i: p2+i;
    if (p == ULONG_MAX) {
      tmpSeq2[i] = 0;
      break;
    }
    tmpSeq2[i] = ((is_rev)? "TGCAN": "ACGTN")[_get_pac(pac, p)];
  }
  for (s = 0; s < maxDepth; ++s) {
    if (tmpSeq1[s] != tmpSeq2[s]) {
      tmpSeq1[s] = 0;
      break;
    }
  }
  for (size_t i = 0; i < s; ++i) {
    tmpSeq2[i] = tmpSeq1[s-i-1];
  }
  tmpSeq2[s] = 0;
  printf("\t\t\t\t%" PRIu64 "-%" PRIu64 " (%p): %s (%c/%c -> %p)\n", getStateInterval(state)->k, getStateInterval(state)->l, state, tmpSeq2, CIGAR[(state->trace & BACKTRACE_MASK) >> BACKTRACE_OFFSET], "ACGT"[state->trace & NUCLEOTIDE_MASK], state->previousState);
}

/******* Smith-Waterman cell type *******/
/**
 * Contains value for backtrace
 */
typedef struct {
  unsigned int value: MAX_SW_COST_SIZE;
  unsigned int backTrace: SW_BACKTRACE_SIZE;
  unsigned int nucleotide: NUCLEOTIDES_BITS;
} sw_cell_t;

/******* Smith-Waterman type *******/
/**
 * Perform the Smith-Waterman algorithm on the reads and the genome.
 * The matrix uses only a band.
 *                  <--- read ---->      original     new
 *                 +-+-+-+-+-+-+-+-+
 *  ^            ^ |X|X|2| | | | | |      +-+-+        +-+
 *  g     # dels | +-+-+-+-+-+-+-+-+      |3|1|        |1|
 *  e            | |X|1|N| | | | | |      +-+-+      +-+-+
 *  n              +-+-+-+-+-+-+-+-+      |2|N|      |3|N|
 *  o   0 error -> |0|N|N| | | | | |      +-+-+      +-+-+
 *  m              +-+-+-+-+-+-+-+-+                 |2|
 *  e            | |1|N|N| | | | | |                 +-+
 *  |      # ins | +-+-+-+-+-+-+-+-+   1: deletion
 *  v            v |2|N|N| | | | | |   2: insertion
 *                 +-+-+-+-+-+-+-+-+   3: (mis)match
 *
 *
 */
typedef struct {
  unsigned short  *readSequence;
  unsigned int     readLength;
  unsigned short **genomeSequences;
  unsigned int    *genomeLengths;
  bwtinterval_t    genomeInterval;
  uint64_t        *poss;
  int             *strands;
  unsigned int     nGenomeSequences;
  unsigned int     alignmentSize;
  sw_cell_t      **matrix;
  state_t         *states;
} sw_t;

void createSW (sw_t *sw, size_t depth) {
  size_t nCols = depth + 1;
  size_t nRows = 2 * parameters->maxNErrors + 1;
  sw->genomeSequences = (unsigned short **) malloc(MAX_SW_N_SEQUENCES * sizeof(unsigned short *));
  sw->genomeLengths   = (unsigned int *)    calloc(MAX_SW_N_SEQUENCES,  sizeof(unsigned int));
  sw->poss            = (uint64_t *)        malloc(MAX_SW_N_SEQUENCES * sizeof(uint64_t));
  sw->strands         = (int *)             malloc(MAX_SW_N_SEQUENCES * sizeof(int));
  sw->readSequence    = (unsigned short *)  malloc(depth * sizeof(unsigned short));
  for (unsigned int i = 0; i < MAX_SW_N_SEQUENCES; ++i) {
    sw->genomeSequences[i] = (unsigned short *) malloc(depth * sizeof(unsigned short));
  }
  sw->matrix = (sw_cell_t **) malloc(nCols * sizeof(unsigned int *));
  for (size_t i = 0; i < nCols; ++i) {
    sw->matrix[i] = (sw_cell_t *) malloc(nRows * sizeof(unsigned int));
  }
  sw->matrix[0][parameters->maxNErrors].value = 0;
  for (unsigned int i = 1; i <= parameters->maxNErrors; ++i) {
    sw->matrix[0][i+parameters->maxNErrors].value = i;
    sw->matrix[0][i+parameters->maxNErrors].backTrace = SW_DELETION;
  }
  for (unsigned int i = 1; i <= parameters->maxNErrors; ++i) {
    sw->matrix[i][parameters->maxNErrors-i].value = i;
    sw->matrix[i][parameters->maxNErrors-i].backTrace = SW_INSERTION;
  }
  sw->states = (state_t *) malloc(nCols * sizeof(state_t));
}

void unsetReadSequence (sw_t *sw) {
  sw->readLength = 0;
}

void addReadSequence (sw_t *sw, sequence_t sequence, unsigned int length) {
  //TODO double-check the direction
  //printf("\t\t\tSetting read sequence (%d): ", length);
  //printSequence(sequence, length);
  //printf("\n");
  for (unsigned int i = 0; i < length; ++i) {
    sw->readSequence[sw->readLength] = sequence & NUCLEOTIDE_MASK;
    sequence >>= NUCLEOTIDES_BITS;
    if (sw->readLength < parameters->maxNErrors) {
      sw->matrix[sw->readLength+1][parameters->maxNErrors-sw->readLength-1].nucleotide = sw->readSequence[sw->readLength];
    }
    ++sw->readLength;
  }
  //printf("\t\t\tSetting read sequence (%d): ", length);
  //printSequenceLong(sw->readSequence, sw->readLength);
  //printf("\n");
}

unsigned int setGenomeSequences (sw_t *sw, state_t *state) {
  assert(state->interval.k <= state->interval.l);
  int      isRev;
  bwtint_t pos;
  //unsigned int nSequencesUndup = 0;
  sw->nGenomeSequences = state->interval.l - state->interval.k + 1;
  // copy all sequences in the interval
  for (unsigned int sequenceId = 0; sequenceId < sw->nGenomeSequences; ++sequenceId) {
    pos = bwt_sa(bwt, state->interval.k + sequenceId);
    pos = bns_depos(bns, pos, &isRev);
    //printf("Pos: %lu/%lu/%lu\n", pos, bwt->bwt_size, bwt->seq_len); fflush(stdout);
    //bns_pos2rid(bns, pos);
    //printf("\t\t\tSetting genome sequence: ");
    sw->poss[sequenceId]          = pos;
    sw->strands[sequenceId]       = isRev;
    sw->genomeLengths[sequenceId] = 0;
    for (unsigned int i = 0; i < MAX_SW_LENGTH; ++i) {
      pos = (isRev)? pos+1: pos-1;
      uint64_t c = DNA5_TO_INT_REV[isRev][_get_pac(pac, pos)];
      sw->genomeSequences[sequenceId][i] = c;
      //printf("%c", "ACGT"[c]);
      ++sw->genomeLengths[sequenceId];
      if ((isRev) && (pos == bwt->seq_len)) {
        break;
      }
      else if ((! isRev) && (pos == 0)) {
        break;
      }
    }
    //printf("\n"); fflush(stdout);
    //printf("\t\t\tTransformed into:        ");
    //printSequence(sw->genomeSequences[sequenceId], sw->genomeLength);
    //printf("\n"); fflush(stdout);
  }
  sw->genomeInterval = state->interval;
  //TODO remove duplicate sequences
  /*
  for (unsigned int sequenceId = 1; sequenceId < sw->nGenomeSequences; ++sequenceId) {
    if (memcmp(sw->genomeSequences[sequenceId], sw->genomeSequences[nSequencesUndup], sw->genomeLengths[sequenceId] * sizeof(unsigned short)) == 0) {
      sw->genomeSequences[sequenceId] = NULL;
      ++nSequencesUndup;
    }
  }
  */
  return sw->nGenomeSequences;
}

void setGenomeSequence (sw_t *sw, unsigned int genomeSequenceId) {
  //printf("\t\t\tUsing genome sequence: ");
  //printSequenceLong(sw->genomeSequences[genomeSequenceId], sw->genomeLengths[genomeSequenceId]);
  //printf("\n");
  for (unsigned int i = 1; i <= MIN(parameters->maxNErrors, sw->genomeLengths[genomeSequenceId]); ++i) {
    sw->matrix[0][parameters->maxNErrors+i].nucleotide = getNucleotide(sw->genomeSequences[genomeSequenceId][i], i-1);
    //printf("\t\tinsertion (0, %zu): %c\n", parameters->maxNErrors+i, "ACGT"[sw->matrix[0][parameters->maxNErrors+i].nucleotide]);
  }
}

bool tryNoDiffSW (sw_t *sw, unsigned int genomeSequenceId) {
  //printf("\t\t\tShrinking genome to:     ");
  //printSequence(genomeSequence, sw->readLength);
  //printf("\n");
  unsigned short genomeNucleotide;
  if (sw->readLength != sw->genomeLengths[genomeSequenceId]) {
    return false;
  }
  if (memcmp(sw->readSequence, sw->genomeSequences[genomeSequenceId], sw->readLength * sizeof(unsigned short)) != 0) {
    return false;
  }
  sw->alignmentSize = sw->readLength;
  //printf("\t\tFound right away!\n");
  for (unsigned int i = 0; i < sw->genomeLengths[genomeSequenceId]; ++i) {
    genomeNucleotide = sw->genomeSequences[genomeSequenceId][i];
    sw->states[i].trace = MATCH | genomeNucleotide;
    sw->states[i].previousState = NULL;
    //printf("\t\t\tState @ %u\n", i);
    //printState(&sw->states[i], 100);
    //printf("\n");
  }
  return true;
}

unsigned int getScore (sw_t *sw, unsigned int genomeSequenceId, unsigned int maxNErrors) {
  uint64_t reconstructedGenomeSequence = 0;
  unsigned short readNucleotide, genomeNucleotide;
  unsigned int readId, genomeId, xId, yId, startingYId = 2 * maxNErrors + 2, bestYId = startingYId;
  unsigned int v1, v2, v3;
  unsigned int minValue = maxNErrors + 1;
  unsigned int reversedBacktraceId;
  bwtinterval_t previousInterval;
  bool match;
  if (tryNoDiffSW(sw, genomeSequenceId)) {
    return 0;
  }
  if (maxNErrors == 0) {
    return 1;
  }
  for (unsigned int i = 1; i <= maxNErrors+1; ++i) {
    sw->matrix[0][i].nucleotide = sw->genomeSequences[genomeSequenceId][i-1];
  }
  // local alignment
  //printf("\t\tAt most: %u errors\n", maxNErrors); fflush(stdout);
  for (readId = 0; readId < sw->readLength; ++readId) {
    readNucleotide = sw->readSequence[readId];
    minValue = maxNErrors + 1;
    for (int diff = -maxNErrors; diff <= (int) maxNErrors; ++diff) {
      if (((int) readId) + diff >= 0) {
        genomeId = readId + diff;
        yId = diff + ((int) parameters->maxNErrors);
        if (genomeId < sw->genomeLengths[genomeSequenceId]) {
          genomeNucleotide = sw->genomeSequences[genomeSequenceId][genomeId];
          match = (readNucleotide == genomeNucleotide);
          v1 = (diff == (int) -maxNErrors)? maxNErrors+1: (unsigned int) sw->matrix[readId+1][yId-1].value + 1;
          v2 = (diff == (int) maxNErrors)?  maxNErrors+1: (unsigned int) sw->matrix[readId][yId+1].value + 1;
          v3 = sw->matrix[readId][yId].value + (match? 0: 1);
          //printf("\t\t\t\t\t%u\n", v1); fflush(stdout);
          //printf("\t\t\t\t\t%u\n", v2); fflush(stdout);
          //printf("\t\t\t\t\t%u\n", v3); fflush(stdout);
          //printf("\t\t\t\t(%u, %i/%i) Values %u %u %u / %u\n", readId, diff, yId, v1, v2, v3, maxNErrors); fflush(stdout);
          if ((v3 <= v1) && (v3 <= v2)) {
            sw->matrix[readId+1][yId].value = v3;
            sw->matrix[readId+1][yId].backTrace = (match)? SW_MATCH: SW_MISMATCH;
            sw->matrix[readId+1][yId].nucleotide = genomeNucleotide;
            //printf("\t\t\t\t(mis)match (%u, %u) <- %u (%c)\n", readId+1, yId, v3, "ACGT"[readNucleotide]);
          }
          else if (v1 <= v2) {
            sw->matrix[readId+1][yId].value = v1;
            sw->matrix[readId+1][yId].backTrace = SW_DELETION;
            sw->matrix[readId+1][yId].nucleotide = genomeNucleotide;
            //printf("\t\t\t\tdeletion (%u, %u) <- %u (%c)\n", readId+1, yId, v1, "ACGT"[genomeNucleotide]);
          }
          else {
            sw->matrix[readId+1][yId].value = v2;
            sw->matrix[readId+1][yId].backTrace = SW_INSERTION;
            sw->matrix[readId+1][yId].nucleotide = readNucleotide;
            //printf("\t\t\t\tinsertion (%u, %u) <- %u (%c)\n", readId+1, yId, v2, "ACGT"[readNucleotide]);
          }
          minValue = MIN(minValue, sw->matrix[readId+1][yId].value);
        }
        else {
          sw->matrix[readId+1][yId].value = maxNErrors + 1;
          sw->matrix[readId+1][yId].backTrace = SW_INSERTION;
          sw->matrix[readId+1][yId].nucleotide = readNucleotide;
        }
      }
    }
    //printf("\t\tminValue %u @ %u\n", minValue, readId); fflush(stdout);
    if (minValue == maxNErrors + 1) {
      return maxNErrors + 1;
    }
  }
  //printf("\t\tFound with mismatches\n");
  // get end of backtrace, favoring (mis)matches
  readId = sw->readLength;
  //printf("\t\t\t\tmax #errors %u, readId %u\n", maxNErrors, readId);
  for (unsigned int diff = 0; diff <= maxNErrors; ++diff) {
    //printf("\t\t\t\tdiff %u\n", diff); fflush(stdout);
    if (diff == 0) {
      yId = parameters->maxNErrors;
      //printf("\t\t\t\t\ttest 0: (%u, %u) = %u\n", readId, yId, sw->matrix[readId][yId].value); fflush(stdout);
      assert(sw->matrix[readId][yId].value >= minValue);
      if (sw->matrix[readId][yId].value == minValue) {
        bestYId = yId;
        break;
      }
    }
    else {
      if (readId >= diff) {
        yId = parameters->maxNErrors - diff;
        //printf("\t\t\t\t\ttest -1: (%u, %u) = %u\n", readId, yId, sw->matrix[readId][yId].value); fflush(stdout);
        assert(sw->matrix[readId][yId].value >= minValue);
        if (sw->matrix[readId][yId].value == minValue) {
          bestYId = yId;
          break;
        }
      }
      yId = parameters->maxNErrors + diff;
      //printf("\t\t\t\t\ttest +1: (%u, %u) = %u\n", readId, yId, sw->matrix[readId][yId].value); fflush(stdout);
      assert(sw->matrix[readId][yId].value >= minValue);
      if (sw->matrix[readId][yId].value == minValue) {
        bestYId = yId;
        break;
      }
    }
  }
  //printf("\t\tbest/starting/min: %u/%u/%u\n", bestYId, startingYId, minValue);
  assert(bestYId != startingYId);
  // get size of backtracing
  xId = readId;
  yId = bestYId;
  sw->alignmentSize = 0;
  while((xId != 0) || (yId != parameters->maxNErrors)) {
    ++sw->alignmentSize;
    switch(sw->matrix[xId][yId].backTrace) {
      case SW_MATCH:
      case SW_MISMATCH:
        //printf("\t\t\t\t(mis)match (%i, %i)\n", xId, yId); fflush(stdout);
        assert(xId > 0);
        --xId;
        break;
      case SW_INSERTION:
        //printf("\t\t\t\tinsertion (%i, %i)\n", xId, yId); fflush(stdout);
        assert(xId > 0);
        --xId;
        ++yId;
        break;
      case SW_DELETION:
        //printf("\t\t\t\tdeletion (%i, %i)\n", xId, yId); fflush(stdout);
        assert(yId > 0);
        --yId;
        break;
      default:
        assert(false);
    }
  }
  //printf("\t\t\tdone, size: %u\n", sw->alignmentSize); fflush(stdout);
  // backtrace
  xId = readId;
  yId = bestYId;
  for (unsigned int backtraceId = 0; backtraceId < sw->alignmentSize; ++backtraceId) {
    reversedBacktraceId = sw->alignmentSize - backtraceId - 1;
    sw->states[reversedBacktraceId].trace = (sw->matrix[xId][yId].backTrace << BACKTRACE_OFFSET) | sw->matrix[xId][yId].nucleotide;
    switch(sw->matrix[xId][yId].backTrace) {
      case SW_MATCH:
      case SW_MISMATCH:
        //printf("\t\t\t\t(mis)match (%u, %u) %c\n", xId, yId, "ACGT"[sw->matrix[xId][yId].nucleotide]); fflush(stdout);
        reconstructedGenomeSequence <<= NUCLEOTIDES_BITS;
        reconstructedGenomeSequence |= sw->matrix[xId][yId].nucleotide;
        --xId;
        break;
      case SW_INSERTION:
        //printf("\t\t\t\tinsertion (%u, %u) %c\n", xId, yId, "ACGT"[sw->matrix[xId][yId].nucleotide]); fflush(stdout);
        --xId;
        ++yId;
        break;
      case SW_DELETION:
        //printf("\t\t\t\tdeletion (%u, %u)\n", xId, yId); fflush(stdout);
        reconstructedGenomeSequence <<= NUCLEOTIDES_BITS;
        reconstructedGenomeSequence |= sw->matrix[xId][yId].nucleotide;
        --yId;
        break;
      default:
        assert(false);
    }
  }
  assert(xId == 0);
  assert(yId == parameters->maxNErrors);
  // add BWT intervals
  previousInterval = sw->genomeInterval;
  for (unsigned int backtraceId = 0; backtraceId < sw->alignmentSize; ++backtraceId) {
  //for (int backtraceId = sw->alignmentSize - 1; backtraceId >= 0; --backtraceId) {
    //printf("\t\t\tbacktrace: %u\n", sw->states[backtraceId].trace);
    if ((sw->states[backtraceId].trace & BACKTRACE_MASK) != INSERTION) {
      genomeNucleotide = reconstructedGenomeSequence & NUCLEOTIDE_MASK;
      /*
      bwt_2occ(bwt, previousInterval.k-1, previousInterval.l, genomeNucleotide, &sw->states[backtraceId].interval.k, &sw->states[backtraceId].interval.l);
      sw->states[backtraceId].interval.k = bwt->L2[genomeNucleotide] + sw->states[backtraceId].interval.k + 1;
      sw->states[backtraceId].interval.l = bwt->L2[genomeNucleotide] + sw->states[backtraceId].interval.l;
      */
      sw->states[backtraceId].previousState = NULL;
      reconstructedGenomeSequence >>= NUCLEOTIDES_BITS;
    }
    else {
      sw->states[backtraceId].interval = previousInterval;
      sw->states[backtraceId].previousState = NULL;
    }
    //printf("\t\t\tState @ %u\n", backtraceId);
    //printState(&sw->states[backtraceId], 100);
    //printf("\n"); fflush(stdout);
    previousInterval = sw->states[backtraceId].interval;
  }
  //printf("\t\t\tgenome sequence: %lu\n", reconstructedGenomeSequence); fflush(stdout);
  assert(reconstructedGenomeSequence == 0);
  return minValue;
}


typedef struct {
  unsigned int   index;    
  bwtinterval_t  previousIntervals[BWT_BUFFER_SIZE];
  unsigned short nucleotides[BWT_BUFFER_SIZE];
  bwtinterval_t  nextIntervals[BWT_BUFFER_SIZE];
  bool           full;
} bwt_buffer_t;

void createBwtBuffer (bwt_buffer_t *bwtBuffer) {
  bwtBuffer->index = 0;
  bwtBuffer->full  = false;
}

void addToBwtBuffer (bwt_buffer_t *bwtBuffer, bwtinterval_t previousInterval, unsigned short nucleotide, bwtinterval_t nextInterval) {
  bwtBuffer->previousIntervals[bwtBuffer->index] = previousInterval;
  bwtBuffer->nucleotides[bwtBuffer->index] = nucleotide;
  bwtBuffer->nextIntervals[bwtBuffer->index] = nextInterval;
  ++bwtBuffer->index;
  if (bwtBuffer->index == BWT_BUFFER_SIZE) {
    bwtBuffer->index = 0;
    bwtBuffer->full = true;
  }
}

bwtinterval_t *findInBwtBuffer (bwt_buffer_t *bwtBuffer, bwtinterval_t previousInterval, unsigned short nucleotide) {
  return NULL;
  unsigned int searchIndex = 0;
  ++stats->nBufferCalls;
  if (bwtBuffer->index == 0) {
    if (! bwtBuffer->full) {
      return NULL;
    }
    searchIndex = BWT_BUFFER_SIZE - 1;
  }
  while (searchIndex != bwtBuffer->index) {
    if (compareBwtIntervals(&previousInterval, &bwtBuffer->previousIntervals[searchIndex]) && (nucleotide == bwtBuffer->nucleotides[searchIndex])) {
      ++stats->nBufferCallSucesses;
      return &bwtBuffer->nextIntervals[searchIndex];
    }
    if (searchIndex == 0) {
      if (! bwtBuffer->full) {
        return NULL;
      }
      searchIndex = BWT_BUFFER_SIZE - 1;
    }
    else { 
      --searchIndex;
    }
  }
  return NULL;
}

bool goDownBwt (bwt_buffer_t *bwtBuffer, state_t *previousState, unsigned short nucleotide, bwtinterval_t *newInterval) {
  //printf("    Going down BWT from range %" PRId64 "-%" PRIu64 " and nt %hu, preprocessed: %s\n", getStateInterval(previousState)->k, getStateInterval(previousState)->l, nucleotide, (previousState->trace & PREPROCESSED)? "true": "false");
  bwtinterval_t *newIntervalInBuffer = findInBwtBuffer(bwtBuffer, previousState->interval, nucleotide);
  ++stats->nDown;
  if (newIntervalInBuffer == NULL) {
    bwt_2occ(bwt, previousState->interval.k-1, previousState->interval.l, nucleotide, &newInterval->k, &newInterval->l);
    newInterval->k = bwt->L2[nucleotide] + newInterval->k + 1;
    newInterval->l = bwt->L2[nucleotide] + newInterval->l;
  }
  else {
    *newInterval = *newIntervalInBuffer;
  }
  //addToBwtBuffer(bwtBuffer, previousState->interval, nucleotide, newInterval);
  //if (newState->interval.k <= newState->interval.l) printf("      ok!\n");
  return (newInterval->k <= newInterval->l);
}


/******* States type *******/
/**
 * The states stores the mapping information of the prefixes of reads.
 * It is:
 *   - an array of states: 1st dimension is the depth, 2nd is the number of errors
 *   - the number of states per depth, per errors
 *   - the number of states per depth
 *   - the min number of error per depth
 *   - the max number of error per depth
 *   - the tree size
 */

typedef struct {
  state_t    ***states;
  size_t      **nStates;
  size_t       *nStatesPerPosition;
  size_t       *minErrors;
  size_t       *maxErrors;
  size_t       depth;
  sw_t         *sw;
  bwt_buffer_t *bwtBuffer;
} states_t;

/*
bool goDownBwt (state_t *previousState, unsigned short nucleotide, state_t *newState) {
  ++stats->nDown;
  //printf("    Going down BWT from range %" PRId64 "-%" PRIu64 " and nt %hu, preprocessed: %s\n", getStateInterval(previousState)->k, getStateInterval(previousState)->l, nucleotide, (previousState->trace & PREPROCESSED)? "true": "false");
  bwt_2occ(bwt, getStateInterval(previousState)->k-1, getStateInterval(previousState)->l, nucleotide, &newState->interval.k, &newState->interval.l);
  newState->interval.k = bwt->L2[nucleotide] + newState->interval.k + 1;
  newState->interval.l = bwt->L2[nucleotide] + newState->interval.l;
  //if (newState->interval.k <= newState->interval.l) printf("      ok!\n");
  newState->trace = 0;
  return (newState->interval.k <= newState->interval.l);
}
*/

void printStates (states_t *states, size_t depth) {
  for (size_t i = 0; i < depth; ++i) {
    if (states->minErrors[i] == SIZE_MAX) {
      printf("\t\t\t\t(%zu) ---\n", i);
    }
    else {
      printf("\t\t\t\t(%zu)", i);
      for (size_t j = states->minErrors[i]; j <= states->maxErrors[i]; ++j) {
        printf(" %zu errors: %zu states ", j, states->nStates[i][j]);
      }
      printf("\n");
    }
  }
}

unsigned int getParentStateId (unsigned int stateId) {
  return (stateId-1) / 2;
}

unsigned int getLowerStateId (unsigned int stateId) {
  return ((2 * stateId) + 1);
}

unsigned int getGreaterStateId (unsigned int stateId) {
  return ((2 * stateId) + 2);
}

void swapStates(state_t *s1, state_t *s2) { 
  state_t tmp = *s1; 
  *s1 = *s2; 
  *s2 = tmp; 
} 

bool areStatesEqual (state_t *state1, state_t *state2) {
  return ((getStateInterval(state1)->k == getStateInterval(state2)->k) && (getStateInterval(state1)->l == getStateInterval(state2)->l));
}

int stateComp (const state_t *state1, const state_t *state2) {
  int i = getStateInterval((state_t *) state1)->k - (getStateInterval((state_t *) state2)->k);
  if (i != 0) {
    return i;
  }
  return (getStateInterval((state_t *) state1)->k - (getStateInterval((state_t *) state2)->k));
}

int sortCompareStates (const void *state1, const void *state2) {
  return (getStateInterval((state_t *) state1)->k - (getStateInterval((state_t *) state2)->k));
}

bool canMerge (state_t *state1, state_t *state2) {
  if ((getStateInterval(state1)->l+1 < getStateInterval(state2)->k) || (getStateInterval(state1)->k+1 > getStateInterval(state2)->l)) {
    return false;
  }
  getStateInterval(state1)->k = MIN(getStateInterval(state1)->k, getStateInterval(state2)->k);
  getStateInterval(state1)->l = MAX(getStateInterval(state1)->l, getStateInterval(state2)->l);
  return true;
}

void simplifyStates (states_t *states, size_t depth, size_t nErrors) {
  //return nStates;
  assert(depth < states->depth);
  assert(nErrors <= states->maxErrors[depth]);
  size_t previousNStates = states->nStates[depth][nErrors];
  size_t nextNStates = 0;
  state_t *theseStates = states->states[depth][nErrors];
  if (previousNStates <= 1) {
    return;
  }
  //TODO check that the pointers are good!
  //printf("\t\t\t\tSimplify states from %zu ", nStates);
  //printf("Entering Simplify States @ depth %zu with %zu errors and %zu elements.\n", depth, nErrors, previousNStates); fflush(stdout);
  qsort(theseStates, previousNStates, sizeof(state_t), sortCompareStates);
  for (size_t secondStateId = 1; secondStateId < previousNStates; ++secondStateId) {
    //printf("\tCurrent state: %zu/%zu/%zu\n", nextNStates, secondStateId, previousNStates); fflush(stdout);
    if (! areStatesEqual(&theseStates[nextNStates], &theseStates[secondStateId])) {
    //if (! canMerge(&states[firstStateId], &states[secondStateId])) {
      ++nextNStates;
      if (nextNStates < secondStateId) {
        theseStates[nextNStates] = theseStates[secondStateId];
      }
      assert(nextNStates <= previousNStates);
      assert(nextNStates <= secondStateId);
      assert(secondStateId < N_STATES);
    }
  }
  ++nextNStates;
  //printf("to %zu\n", firstStateId+1);
  states->nStates[depth][nErrors] = nextNStates;
  assert(states->nStatesPerPosition[depth] >= previousNStates - nextNStates);
  states->nStatesPerPosition[depth] -= previousNStates - nextNStates;
}

void heapifyStates (state_t *states, size_t nStates) {
  return;
  size_t currentId = nStates-1;
  size_t parentId  = getParentStateId(currentId);
  while ((currentId != 0) && (stateComp(&states[parentId], &states[currentId]) > 0)) { 
    swapStates(&states[currentId], &states[parentId]); 
    currentId = parentId;
  }
}

bool isStateInserted (state_t *states, size_t nStates, state_t *state) {
  return false;
  size_t currentId = 0;
  int cmp;
  while (currentId < nStates) {
    cmp = stateComp(state, &states[currentId]);
    if (cmp == 0) {
      return true;
    }
    else if (cmp < 0) {
      currentId = getLowerStateId(currentId);
    }
    else {
      currentId = getGreaterStateId(currentId);
    }
  }
  return false;
}

state_t *addState (states_t *states, size_t depth, size_t nErrors) {
  assert(depth <= states->depth);
  //printf("\t\t\tAdding one state (%" PRIu64 ", %" PRIu64 ") at (depth = %zu, # errors = %zu), %zu/%d occupied\n", state->interval.k, state->interval.l, depth, nErrors, states->nStates[depth][nErrors], N_STATES); fflush(stdout);
  //printState(state, states->depth);
  /*
  if (isStateInserted(states->states[depth][nErrors], states->nStates[depth][nErrors], state)) {
    return NULL;
  }
  */
  if (states->nStates[depth][nErrors] == N_STATES-1) {
    printf("Exiting because # states = %zu >= %i at depth %zu with %zu errors, before simplification.\n", states->nStates[depth][nErrors], N_STATES, depth, nErrors);
    printStates(states, depth);
    states->nStates[depth][nErrors] = N_STATES;
    return NULL;
  }
  //states->states[depth][nErrors][states->nStates[depth][nErrors]] = *state;
  ++states->nStates[depth][nErrors];
  ++states->nStatesPerPosition[depth];
  stats->maxNStates = MAX(stats->maxNStates, states->nStates[depth][nErrors]);
  heapifyStates(states->states[depth][nErrors], states->nStates[depth][nErrors]);
  if (states->minErrors[depth] == SIZE_MAX) {
    states->minErrors[depth] = nErrors;
    states->maxErrors[depth] = nErrors;
  }
  else {
    if (nErrors < states->minErrors[depth]) {
      states->minErrors[depth] = nErrors;
    }
    if (nErrors > states->maxErrors[depth]) {
      states->maxErrors[depth] = nErrors;
    }
  }
  return &states->states[depth][nErrors][states->nStates[depth][nErrors]-1];
}

states_t *initializeStates(size_t treeSize) {
  states_t *states           = (states_t *)     malloc(sizeof(states_t));
  states->depth              = treeSize + 1 + parameters->maxNErrors;
  states->states             = (state_t ***)    malloc(states->depth * sizeof(states_t **));
  states->nStates            = (size_t **)      malloc(states->depth * sizeof(size_t *));
  states->nStatesPerPosition = (size_t *)       calloc(states->depth,  sizeof(size_t));
  states->minErrors          = (size_t *)       malloc(states->depth * sizeof(size_t));
  states->maxErrors          = (size_t *)       malloc(states->depth * sizeof(size_t));
  states->sw                 = (sw_t *)         malloc(sizeof(sw_t));
  states->bwtBuffer          = (bwt_buffer_t *) malloc(sizeof(bwt_buffer_t));
  for (size_t depth = 0; depth < states->depth; ++depth) {
    states->states[depth]    = (state_t **) malloc((parameters->maxNErrors+1) * sizeof(states_t *));
    states->nStates[depth]   = (size_t *)   calloc((parameters->maxNErrors+1),  sizeof(size_t));
    states->minErrors[depth] = SIZE_MAX;
    states->maxErrors[depth] = SIZE_MAX;
    for (size_t nErrors = 0; nErrors <= parameters->maxNErrors; ++nErrors) {
      states->states[depth][nErrors] = (state_t *) malloc(N_STATES * sizeof(states_t));
    }
  }
  state_t *baseState = addState(states, 0, 0);
  baseState->interval.k    = 0;
  baseState->interval.l    = bwt->seq_len;
  baseState->trace         = 0;
  baseState->previousState = NULL;
  createSW(states->sw, states->depth);
  createBwtBuffer(states->bwtBuffer);
  return states;
}

void backtrackStates(states_t *states, size_t level) {
  //printf("\t\t\tBacktracking to level %zu\n", level);
  for (size_t i = level; i <= states->depth; ++i) {
    if (states->nStatesPerPosition[i] == 0) {
      return;
    }
    states->nStatesPerPosition[i] = 0;
    states->minErrors[i]          = SIZE_MAX;
    states->maxErrors[i]          = SIZE_MAX;
    for (size_t j = 0; j <= parameters->maxNErrors; ++j) {
      states->nStates[i][j] = 0;
    }
  }
}

void freeStates(states_t *states) {
  for (size_t i = 0; i < states->depth; ++i) {
    for (size_t j = 0; j <= parameters->maxNErrors; ++j) {
      free(states->states[i][j]);
    }
    free(states->states[i]);
  }
  free(states->states);
  free(states->nStates);
  free(states->nStatesPerPosition);
  free(states->minErrors);
  free(states->maxErrors);
  free(states);
}

typedef struct {
  bool         isSet;
  size_t       depth;
  int          isRev;
  int          rid;
  unsigned int nMisses;
  bwtint_t     pos;
} shortCut_t;

void unsetShortCut (shortCut_t *shortCut) {
  shortCut->depth   = DEPTH_SHORT_CUT;
  shortCut->isSet   = false;
  shortCut->nMisses = 0;
}

shortCut_t *initializeShortCut () {
  shortCut_t *shortCut = (shortCut_t *) malloc(sizeof(shortCut_t));
  unsetShortCut(shortCut);
  return shortCut;
}

void setShortCut (shortCut_t *shortCut, bwtint_t k, size_t depth) {
  shortCut->isSet = true;
  shortCut->depth = depth;
  shortCut->pos   = bwt_sa(bwt, k);
  shortCut->pos   = bns_depos(bns, shortCut->pos, &shortCut->isRev);
  shortCut->rid   = bns_pos2rid(bns, shortCut->pos);
}

void resetShortCut (shortCut_t *shortCut, size_t depth) {
  shortCut->pos = (shortCut->isRev)? shortCut->pos+(depth-shortCut->depth): shortCut->pos-(depth-shortCut->depth);
  shortCut->depth = depth;
}

void incShortCut (shortCut_t *shortCut) {
  shortCut->pos = (shortCut->isRev)? shortCut->pos+1: shortCut->pos-1;
  ++shortCut->depth;
}

void addMiss (shortCut_t *shortCut) {
  shortCut->depth += shortCut->nMisses;
  ++shortCut->nMisses;
}

/******* Path type *******/
/**
 * The path stores the path to a cell of the reads tree.
 * It is:
 *   - the nucleotide chose up to that cell
 *   - the corresponding cellId for each cell of the path
 *   - the tree depth
 *   - the corresponding read
 *   - the index of the next char to be written for the read (not the one you want to change when you change the branch)
 *   - nCells is actually the index of the last used cellId
 */

typedef struct {
  unsigned short *nucleotides;
  uint64_t       *cellIds;
  size_t          nCells;
  size_t          maxDepth;
  size_t          depth;
  edge_t         *edges;
  size_t          edgeLength;
  char           *read;
  size_t          readPos;
  shortCut_t     *shortCut;
} path_t;

void printPath (path_t *path) {
  assert(path->depth <= path->maxDepth);
  printf("Read: %s\n", path->read + path->readPos);
  for (size_t i = 0; i < path->depth; ++i) {
    assert(path->nucleotides[i] < N_NUCLEOTIDES);
    printf("%c ", "ACGT"[(int) path->nucleotides[i]]);
  }
  printf("\n");
  for (size_t i = 0; i <= path->nCells; ++i) {
    printf("%" PRIu64 " ", path->cellIds[i]);
    if ((i >= TREE_BASE_SIZE) && (i < path->nCells)) {
      printEdge(&path->edges[i]);
    }
    printf("  ");
  }
  printf("\n(depth: %zu, # cells: %zu, read pos: %zu, cellId: %" PRIu64 ", edge length: %zu)\n", path->depth, path->nCells, path->readPos, path->cellIds[path->nCells], path->edgeLength); fflush(stdout);
}

path_t *initializePath (size_t maxDepth) {
  path_t *path         = (path_t *)         malloc(sizeof(path_t));
  path->nucleotides    = (unsigned short *) malloc(maxDepth * sizeof(unsigned short));
  path->cellIds        = (uint64_t *)       malloc((maxDepth+1) * sizeof(uint64_t));
  path->edges          = (edge_t *)         malloc((maxDepth+1) * sizeof(edge_t));
  path->nCells         = 0;
  path->read           = (char *)           malloc((maxDepth+1) * sizeof(char));
  path->maxDepth       = maxDepth;
  path->depth          = 0;
  path->edgeLength     = 0;
  path->cellIds[0]     = 0;
  path->read[maxDepth] = 0;
  path->readPos        = maxDepth;
  path->shortCut       = initializeShortCut();
  return path;
}

void freePath (path_t *path) {
  free(path->nucleotides);
  free(path->cellIds);
  free(path->shortCut);
  free(path);
}

void appendNucleotidePath (path_t *path, short nucleotide, char c) {
  assert(path->depth < path->maxDepth);
  path->nucleotides[path->depth] = nucleotide;
  ++path->depth;
  --path->readPos;
  path->read[path->readPos] = c;
}

/**
 * Step into the last cells of the tree in a DFS fashion.
 */
bool goDownTreeBase (path_t *path) {
  //printf("\t\tGoing to base\n");
  appendNucleotidePath(path, 0, 'A');
  path->cellIds[path->nCells+1] = path->cellIds[path->nCells] << NUCLEOTIDES_BITS;
  ++path->nCells;
  return true;
}

/**
 * Step into the first cells of the tree in a DFS fashion.
 */
bool goDownTreeNotBase (const tree_t *tree, path_t *path) {
  //printf("\t\tGoing to tree\n");
  edge_t *edge = NULL;
  unsigned short nucleotide;
  //printPath(path);
  assert(path->depth <= tree->depth);
  assert(path->nCells <= tree->depth);
  if (path->edgeLength == 0) {
    cell_t *cell = &tree->cells[path->cellIds[path->nCells]];
    for (nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
      edge = &cell->edges[nucleotide];
      if (isSetEdge(edge)) {
        path->edges[path->nCells] = *edge;
        path->edgeLength          = 1;
        appendNucleotidePath(path, nucleotide, DNA5_TO_CHAR[nucleotide]);
        if (path->edgeLength == edge->length) {
          path->edgeLength = 0;
          path->cellIds[++path->nCells] = edge->cellId;
        }
        //printf("to %zu with read '%s'\n", path->depth, path->read+path->readPos);
        return true;
      }
    }
    //printf("\t\t\tNothing found\n");
    return false;
  }
  edge = &path->edges[path->nCells];
  nucleotide = getEdgeNucleotide(edge, path->edgeLength);
  ++path->edgeLength;
  appendNucleotidePath(path, nucleotide, DNA5_TO_CHAR[nucleotide]);
  if (path->edgeLength == edge->length) {
    path->edgeLength = 0;
    path->cellIds[++path->nCells] = edge->cellId;
  }
  return true;
}

/**
 * Step into the tree in a DFS fashion.
 * Return false if search is exhausted.
 */
bool goDownTree (const tree_t *tree, path_t *path) {
  //printf("\tGo down from height %zu with read '%s'\n", path->depth, path->read+path->readPos);
  if (path->depth < TREE_BASE_SIZE) {
    return goDownTreeBase(path);
  }
  return goDownTreeNotBase(tree, path);
}

/**
 * Step into the last cells of the tree in a BFS fashion.
 */
bool goRightTreeBase (path_t *path) {
  //uint64_t cellId = path->cellIds[path->nCells], nextCellId = cellId + 1, mask = NUCLEOTIDE_MASK;
  unsigned short newNucleotide = 0;
  assert(path->depth <= TREE_BASE_SIZE);
  assert(path->nCells <= TREE_BASE_SIZE);
  //printf("    Entering go right base will cell %" PRIu64 " at depth %zu, last nt %i, and read pos %zu\n", path->cellIds[path->nCells], path->depth, path->nucleotides[path->depth-1], path->readPos);
  for (; path->nucleotides[path->depth-1] == N_NUCLEOTIDES - 1; --path->depth, ++path->readPos) {
    if (path->depth == 1) return false;
  }
  newNucleotide = path->nucleotides[path->depth-1] + 1;
  path->read[path->readPos] = DNA5_TO_CHAR[newNucleotide];
  path->nucleotides[path->depth-1] = newNucleotide;
  /*
  for (unsigned int offset = 0; (cellId & mask) != (nextCellId & mask); ++offset, mask <<= NUCLEOTIDES_BITS) {
    newNucleotide = ((nextCellId & mask) >> (NUCLEOTIDES_BITS * offset));
    //printf("      Offset: %u, mask: %" PRIu64 ", depth: %zu, new char: %c\n", offset, mask, path->depth, DNA5_TO_CHAR[newNucleotide]);
    path->read[path->readPos+offset] = DNA5_TO_CHAR[newNucleotide];
    path->nucleotides[path->depth-offset-1] = newNucleotide;
  }
  */
  path->nCells = path->depth;
  path->cellIds[path->nCells] = path->cellIds[path->nCells] + 1;
  //printf("      Leaving base will cell %" PRIu64 ", nucleotide %c, read %s\n", path->cellIds[path->nCells], "ACGT"[newNucleotide], path->read + path->readPos);
  path->edgeLength = 0;
  return true;
}

/**
 * Step into the first cells of the tree in a BFS fashion.
 */
bool goRightTreeNotBase (const tree_t *tree, path_t *path) {
  assert(path->depth <= path->maxDepth);
  assert(path->nCells <= path->maxDepth);
  edge_t *edge;
  //printf("    not base\n");
  while (true) {
    if (path->depth <= TREE_BASE_SIZE) {
      ++path->nCells;
      return goRightTreeBase(path);
    }
    //printf("    ... trying depth %zu, # cell %zu, edge len %zu\n", path->depth, path->nCells, path->edgeLength);
    assert(path->depth >= path->edgeLength);
    if (path->edgeLength == 0) {
      --path->nCells;
      path->edgeLength = path->edges[path->nCells].length;
    }
    path->depth   -= path->edgeLength;
    path->readPos += path->edgeLength;
    //printf("    ... trying depth %zu, # cell %zu, edge len %zu\n", path->depth, path->nCells, path->edgeLength);
    //printf("    ... cellId: %" PRIu64 ", nt: %i\n", path->cellIds[path->nCells], path->nucleotides[path->depth]);
    //printf("    ");
    //printEdge(&tree->cells[path->cellIds[path->nCells]].edges[path->nucleotides[path->depth]]);
    assert(path->depth < path->maxDepth);
    // following assert does not work if mutation at first nucleotide...
    //assert(tree->cells[path->cellIds[path->nCells]].edges[path->nucleotides[path->depth]].cellId == path->edges[path->nCells].cellId);
    for (size_t nucleotide = path->nucleotides[path->depth]+1; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
      //printf("    ... trying nucleotide %c\n", "ACGTN"[nucleotide]);
      edge = &tree->cells[path->cellIds[path->nCells]].edges[nucleotide];
      //printf("    ... edge ");
      //printEdge(edge);
      //printf("\n");
      if (isSetEdge(edge)) {
        path->nucleotides[path->depth] = nucleotide;
        ++path->depth;
        --path->readPos;
        path->read[path->readPos] = DNA5_TO_CHAR[nucleotide];
        path->edges[path->nCells] = *edge;
        path->edgeLength = 1;
        if (path->edgeLength == edge->length) {
          path->edgeLength = 0;
          path->cellIds[++path->nCells] = edge->cellId;
        }
        return true;
      }
    }
    --path->nCells;
    path->edgeLength = path->edges[path->nCells].length;
  }
  //printf("    ... going to nothing\n");
  return false;
}

/**
 * Step into the tree in a BFS fashion.
 * Return false if search is exhausted.
 */
bool goRightTree (const tree_t *tree, path_t *path) {
  //printf("  Go right read from depth %zu with read '%s'\n", path->depth, path->read+path->readPos);
  unsetShortCut(path->shortCut);
  if (path->depth > TREE_BASE_SIZE) {
    return goRightTreeNotBase(tree, path);
  }
  return goRightTreeBase(path);
}

/**
 * Go to next child node in reads tree.  If none, go to sibling
 */
bool goNextTree (const tree_t *tree, states_t *states, path_t *path, bool mappable) {
  //printf("  Goto next\n");
  if (mappable) {
    //printf("    Mappable\n");
    if (goDownTree(tree, path)) {
      return true;
    }
  }
  //printf("  Going right\n");
  //printStates(states, path->depth);
  //printPath(path);
  if (! goRightTree(tree, path)) {
    return false;
  }
  //printf("  Now\n");
  //printStates(states, path->depth);
  backtrackStates(states, path->depth);
  return true;
}

void writeQname (outputSam_t *outputSam, count_t *counts) {
  size_t qnameLength = sprintf(outputSam->qname, "read%lu_x", ++nReads);
  for (size_t readsFileId = 0; readsFileId < parameters->nReadsFiles; ++readsFileId) {
    qnameLength += sprintf(outputSam->qname+qnameLength, "_%lu", counts[readsFileId]);
  }
}

unsigned int computeMapq (unsigned int nHits, unsigned int nErrors) {
  if ((nHits > 1) || (nErrors >= 40)) {
    return 0;
  }
  else {
    return 40 - nErrors;
  }
}

/**
 * Print one line in the SAM file
 */
void printReadLine (bool forward, unsigned int flag, char *chrName, int64_t pos, bwtint_t nHits, unsigned int hitId, unsigned int nErrors, outputSam_t *outputSam) {
  char *cigar, *seq, *qual;
  if (! outputSam->isBackwardSet) {
    computeReverseComplement(outputSam);
  }
  cigar = (forward)? outputSam->forwardCigar: outputSam->backwardCigar;
  seq   = (forward)? outputSam->forwardSeq:   outputSam->backwardSeq;
  qual  = (forward)? outputSam->forwardQual:  outputSam->backwardQual;
  fprintf(outputSam->file, "%s\t%u\t%s\t%" PRId64 "\t%d\t%s\t*\t0\t0\t%s\t%s\tNH:i:%lu\tHI:i:%u\tIH:i:%lu\tNM:i:%u\n", outputSam->qname, flag, chrName, pos, computeMapq(nHits, nErrors), cigar, seq, qual, nHits, hitId, nHits, nErrors);
}

/**
 * Print a read, which maps only once, and with no error
 */
void printReadUniqueNoError (int strand, char *chrName, int64_t pos, size_t readLength, outputSam_t *outputSam) {
  unsigned int flag = 0;
  setCigarNoError(outputSam, readLength);
  if (strand == 0) {
    flag = CIGAR_REVERSE;
  }
  printReadLine(strand != 0, flag, chrName, pos, 1, 1, 0, outputSam);
  //fprintf(outputSamFile, "%s\t%u\t%s\t%" PRId64 "\t40\t%zuM\t*\t0\t0\t%s\t%s\tNH:i:1\tHI:i:1\tIH:i:1\tNM:i:0\n", qname, flag, chrName, pos, readLength, seq, qual);
}

/**
 * Collect information before printing read
 */
void preparePrintRead (uint64_t pos, int rid, int strand, bwtint_t nHits, unsigned int hitId, unsigned int nErrors, outputSam_t *outputSam) {
  unsigned int flag = (hitId == 0)? 0: CIGAR_SECONDARY_HIT;
  if (strand == 0) {
    flag |= CIGAR_REVERSE;
  }
  printReadLine(strand != 0, flag, bns->anns[rid].name, pos, nHits, hitId, nErrors, outputSam);
}

void computeBacktrace (state_t *state, outputSam_t *outputSam) {
  assert(state != NULL);
  outputSam->backtraceSize = 0;
  while (state->previousState != NULL) {
    //printState(state, 101); fflush(stdout);
    char cigar = CIGAR[state->trace >> BACKTRACE_OFFSET];
    if ((outputSam->backtraceSize == 0) || (outputSam->backtraceCigar[outputSam->backtraceSize-1] != cigar)) {
      outputSam->backtraceCigar[outputSam->backtraceSize] = cigar;
      outputSam->backtraceLengths[outputSam->backtraceSize] = 1;
      ++outputSam->backtraceSize;
    }
    else {
      ++outputSam->backtraceLengths[outputSam->backtraceSize-1];
    }
    state = state->previousState;
  }
}

/**
 * Print all the reads in a given state
 */
void printReadState (state_t *state, size_t depth, bwtint_t nHits, unsigned int nErrors, outputSam_t *outputSam) {
  unsigned int hitId = 0;
  int64_t pos;
  int strand, rid;
  if (nErrors != 0) {
    //printf("depth: %zu\n", depth);
    computeBacktrace(state, outputSam);
    computeCigar(outputSam);
  }
  for (bwtint_t bwtint = getStateInterval(state)->k; bwtint <= getStateInterval(state)->l; ++bwtint) {
    ++hitId;
    pos = bwa_sa2pos(bns, bwt, bwtint, depth, &strand);
    rid = bns_pos2rid(bns, pos);
    pos = pos - bns->anns[rid].offset + 1;
    preparePrintRead(pos, rid, strand, nHits, hitId, nErrors, outputSam);
  }
}

/**
 * Print read in SAM format
 */
void printRead (states_t *states, path_t *path, char *quality, count_t *counts, outputSam_t *outputSam) {
  assert(path->depth <= path->maxDepth);
  size_t depth = path->depth;
  int strand, rid;
  int64_t pos;
  char *seq = path->read + path->readPos;
  char *qual = quality;
  bwtint_t nHits = 0;
  unsigned int nErrors = states->minErrors[depth];
  state_t *theseStates = states->states[depth][nErrors];
  size_t nStates;
  simplifyStates(states, depth, nErrors);
  nStates = states->nStates[depth][nErrors];
  //printf("Quality size: %zu vs %zu\n", strlen(quality), depth);
  //printf("Quality: %s (%p)\n", quality, quality);
  //printf("Read:    %s\n", forwardSeq);
  assert(strlen(quality) == depth);
  //printf("Printing read\n");
  //printPath(path); fflush(stdout);
  //printStates(states, path->depth+1); fflush(stdout);
  writeQname(outputSam, counts);
  outputSam->isBackwardSet = false;
  //printf("depth: %zu/%zu, # states: %zu, nErrors: %u, Seq: %s, qual: %s\n", depth, path->maxDepth, nStates, nErrors, seq, qual); fflush(stdout);
  //printState(&theseStates[nStates-1], path->depth); fflush(stdout);
  memcpy(outputSam->forwardSeq,  seq,  (depth+1) * sizeof(char));
  memcpy(outputSam->forwardQual, qual, (depth+1) * sizeof(char));
  for (size_t i = 0; i < nStates; ++i) {
    nHits += theseStates[i].interval.l - theseStates[i].interval.k + 1;
  }
  if (nHits > parameters->maxNHits) {
    fprintf(outputSam->file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tNH:i:%lu\tNM:i:%u\n", outputSam->qname, seq, qual, nHits, nErrors);
    return;
  }
  if ((nHits == 1) && (nErrors == 0)) {
    pos = bwa_sa2pos(bns, bwt, theseStates[0].interval.k, depth, &strand);
    rid = bns_pos2rid(bns, pos);
    pos = pos - bns->anns[rid].offset + 1;
    printReadUniqueNoError(strand, bns->anns[rid].name, pos, depth, outputSam);
    return;
  }
  //printf("Read: %s, read length: %zu\n", forwardSeq, readLength);
  //printStates(states, depth);
  for (size_t stateId = 0; stateId < nStates; ++stateId) {
    printReadState(&theseStates[stateId], path->depth, nHits, nErrors, outputSam);
  }
}

/**
 * Map without error
 */
bool mapWithoutError (states_t *states, size_t depth, unsigned short nt, size_t nErrors) {
  assert(depth > 0);
  state_t *previousState;
  state_t *nextState;
  bwtinterval_t nextInterval;
  bool mapFound = false;
  //printf("    Mapping %c without error at depth %zu with %zu errors and %zu states\n", "ACGT"[nt], depth, nErrors, states->nStates[depth-1][nErrors]);
  /*
  if (states->nStates[depth-1][nErrors] >= MANY_STATES) {
    states->nStates[depth-1][nErrors] = simplifyStates(states->states[depth-1][nErrors], states->nStates[depth-1][nErrors]);
  }
  */
  for (size_t stateId = 0; stateId < states->nStates[depth-1][nErrors]; ++stateId) {
    previousState = &states->states[depth-1][nErrors][stateId];
    if (goDownBwt(states->bwtBuffer, previousState, nt, &nextInterval)) {
      mapFound = true;
      //printState(&nextState, depth);
      nextState = addState(states, depth, nErrors);
      if (nextState == NULL) {
        //printf("      cannot add state\n");
        return false;
      }
      setState(nextState, &nextInterval, MATCH, previousState);
    }
  }
  //printf("    found map: %s\n", mapFound ? "true" : "false");
  return mapFound;
}

/**
 * Find the mappings with nErrors at depth.
 * Supposes that mappings at depth-1 with nErrors are computed.
 * Supposes that mappings at depth with nErrors-1 are computed.
 */
bool _addError (states_t *states, path_t *path, size_t nErrors, size_t depth) {
  //printf("    depth %zu, %zu errors, nucleotide %c, %zu states, min errors: %zu\n", depth, nErrors, "ACGT"[path->nucleotides[depth-1]],  states->nStates[depth-1][nErrors-1], states->minErrors[depth-1]);
  //printf("      Path is: ");
  //for (size_t i = 0; i < depth; ++i) putchar("ACGT"[path->nucleotides[i]]);
  //putchar('\n');
  //TODO may be skipped?
  bwtinterval_t nextInterval;
  state_t *previousState;
  state_t *nextState;
  if ((states->maxErrors[depth] != SIZE_MAX) && (states->maxErrors[depth] >= nErrors)) {
    //printf("      first case: %zu/%zu/%zu %zu/%i\n", states->maxErrors[depth], nErrors, SIZE_MAX, states->nStates[depth][nErrors], N_STATES);
    return (states->nStates[depth][nErrors] < N_STATES);
  }
  for (size_t stateId = 0; stateId < states->nStates[depth-1][nErrors-1]; ++stateId) {
    previousState = &states->states[depth-1][nErrors-1][stateId];
    //printState(state, path->maxDepth);
    // add insertion
    if (! hasTrace(previousState, DELETION)) {
      nextState = addState(states, depth, nErrors);
      if (nextState == NULL) {
        return false;
      }
      setState(nextState, &previousState->interval, INSERTION, previousState);
    }
    // add mismatches
    for (unsigned short nt = 0; nt < N_NUCLEOTIDES; ++nt) {
      if (goDownBwt(states->bwtBuffer, previousState, nt, &nextInterval)) {
        //addState(states, depth-1, nErrors, &newState);
        if (nt != path->nucleotides[depth-1]) {
          //printState(newState, path->maxDepth);
          nextState = addState(states, depth, nErrors);
          if (nextState == NULL) {
            //printf("      third case\n");
            return false;
          }
          setState(nextState, &nextInterval, MISMATCH | nt, previousState);
        }
      }
    }
  }
  // add deletions
  for (size_t stateId = 0; stateId < states->nStates[depth][nErrors-1]; ++stateId) {
    previousState = &states->states[depth][nErrors-1][stateId];
    if (! hasTrace(previousState, INSERTION)) {
      for (unsigned short nt = 0; nt < N_NUCLEOTIDES; ++nt) {
        if (goDownBwt(states->bwtBuffer, previousState, nt, &nextInterval)) {
          nextState = addState(states, depth, nErrors);
          if (nextState == NULL) {
            return false;
          }
          setState(nextState, &nextInterval, DELETION | nt, previousState);
        }
      }
    }
  }
  if (! mapWithoutError(states, depth, path->nucleotides[depth-1], nErrors)) {
    if (states->nStates[depth][nErrors] == N_STATES-1) {
      //printf("      fourth case\n");
      return false;
    }
  }
  //TODO adapth this
  /*
  if (states->nStates[depth][nErrors] >= MANY_STATES) {
    states->nStates[depth][nErrors] = simplifyStates(states->states[depth][nErrors], states->nStates[depth][nErrors]);
  }
  */
  //states->nStates[depth][nErrors] = simplifyStates(states->states[depth][nErrors], states->nStates[depth][nErrors]);
  if (states->nStates[depth][nErrors] >= N_STATES) {
    //printf("      fifth case: # states = %zu >= %i after simplification.\n", states->nStates[depth][nErrors], N_STATES);
    return false;
  }
  return true;
}

/**
 * Map the last nucleotide, with one more error, at depth path->depth-1.
 */
bool addError (states_t *states, path_t *path) {
  size_t nErrors = states->minErrors[path->depth-1] + 1;
  size_t firstDepth;
  bool firstDepthFound = false;
  // Find the first place where the mappings with nErrors are computed
  for (firstDepth = path->depth; (firstDepth > 1) && (! firstDepthFound); --firstDepth) {
    for (size_t nucleotide = 0; (nucleotide < N_NUCLEOTIDES) && (! firstDepthFound); ++nucleotide) {
      firstDepthFound = (states->nStates[firstDepth][nErrors] > 0);
      //if (firstDepthFound) printf("Found first depth at depth = %zu, # errors = %zu, nucleotide = %zu\n", firstDepth, nErrors, nucleotide);
    }
  }
  if (firstDepthFound) firstDepth += 2;
  //printf("  adding error with %zu @ %zu from %zu\n", nErrors, path->depth, firstDepth);
  for (size_t depth = firstDepth; depth <= path->depth; ++depth) {
    if (! _addError(states, path, nErrors, depth)) {
      //printf("    too many errors\n");
      return false;
    }
  }
  return true;
}

/**
 * Compute the mapping with all possible errors
 */
void mapWithErrors (states_t *states, path_t *path) {
  assert(path->depth <= TREE_BASE_SIZE);
  assert(path->depth > 0);
  mapWithoutError(states, path->depth, path->nucleotides[path->depth-1], states->minErrors[path->depth-1]);
  //TODO: Optimize this
  for (unsigned int nErrors = 1; nErrors <= parameters->maxNErrors; ++nErrors) {
    _addError(states, path, nErrors, path->depth);
  }
  //printf("Mapping with errors\n");
  //printStates(states, path->depth-1);
  //printPath(path);
  //TODO check this
  if (path->depth == TREE_BASE_SIZE) {
    for (size_t nErrors = states->minErrors[TREE_BASE_SIZE]; nErrors <= states->maxErrors[TREE_BASE_SIZE]; ++nErrors) {
      simplifyStates(states, TREE_BASE_SIZE, nErrors);
    }
  }
}

bool shortCutCondition (const states_t *states, const tree_t *tree, const path_t *path) {
  if (path->depth < TREE_BASE_SIZE) return false;
  if (path->edgeLength != 0) return false;
  if (isTerminal(&tree->cells[path->cellIds[path->nCells]])) return false;
  if (! isCellUnbranched(tree, &tree->cells[path->cellIds[path->nCells]])) return false;
  //if (getNChildren(&tree->cells[path->cellIds[path->nCells]]) != 1) return false;
  unsigned int nStates = 0;
  for (unsigned int nErrors = states->minErrors[path->depth]; nErrors <= states->maxErrors[path->depth]; ++nErrors) {
    nStates += states->nStates[path->depth][nErrors];
    if (nStates > MAX_SW_N_STATES) return false;
    for (unsigned int stateId = 0; stateId < states->nStates[path->depth][nErrors]; ++stateId) {
      if (states->states[path->depth][nErrors][stateId].interval.l - states->states[path->depth][nErrors][stateId].interval.k >= MAX_SW_N_SEQUENCES) return false;
    }
  }
  return true;
}

bool tryShortCuts (const tree_t *tree, states_t *states, path_t *path, outputSam_t *outputSam) {
  //printf("    Entering tryShortCuts will cell %" PRIu64 " at depth %zu, last nt %i, and read pos %zu, edge len: %zu\n", path->cellIds[path->nCells], path->depth, path->nucleotides[path->depth-1], path->readPos, path->edgeLength);
  //printPath(path);
  state_t bestStates[2 * MAX_EDGE_LENGTH];
  unsigned int currentNErrors, bestNErrors = parameters->maxNErrors + 1, baseNErrors = 0;
  unsigned int bestStateId = 0;
  unsigned int alignmentSize = 0;
  unsigned int genomeAlignmentSize;
  unsigned int nGenomeSequences;
  unsigned int previousDepth = path->depth;
  unsigned int previousNCells = path->nCells;
  unsigned int previousReadPos = path->readPos;
  uint64_t     bestPos = 0;
  int          bestStrand = 0;
  int          rid;
  char         cigar;
  state_t *previousState;
  state_t *nextState;
  uint64_t cellId = path->cellIds[path->nCells];
  cell_t *cell = &tree->cells[cellId];
  edge_t *edge = getFirstEdge(cell);
  assert(edge != NULL);
  unsetReadSequence(states->sw);
  do {
    assert(edge->length != 0);
    path->edges[path->nCells] = *edge;
    path->cellIds[++path->nCells] = edge->cellId;
    addReadSequence(states->sw, edge->sequence, edge->length);
    cellId = edge->cellId;
    cell = &tree->cells[cellId];
    edge = getFirstEdge(cell);
  }
  while (edge != NULL);
  //TODO store all the best paths instead
  for (unsigned int nErrors = states->minErrors[previousDepth]; (nErrors <= states->maxErrors[previousDepth]) && (nErrors <= bestNErrors); ++nErrors) {
    //printf("      Trying with %u errors\n", nErrors);
    simplifyStates(states, previousDepth, nErrors);
    for (unsigned int stateId = 0; stateId < states->nStates[previousDepth][nErrors]; ++stateId) {
      //printf("      Trying with state #%u/%zu %p\n", stateId, states->nStates[previousDepth][nErrors], &states->nStates[previousDepth][nErrors]);
      //TODO: Do not store all the genome sequences
      nGenomeSequences = setGenomeSequences(states->sw, &states->states[previousDepth][nErrors][stateId]);
      for (unsigned int genomeSequenceId = 0; genomeSequenceId < nGenomeSequences; ++genomeSequenceId) {
        setGenomeSequence(states->sw, genomeSequenceId);
        currentNErrors = getScore(states->sw, genomeSequenceId, bestNErrors - nErrors - 1);
        if (currentNErrors + nErrors < bestNErrors) {
          bestNErrors   = currentNErrors + nErrors;
          baseNErrors   = nErrors;
          bestStateId   = stateId;
          alignmentSize = states->sw->alignmentSize;
          bestPos       = states->sw->poss[genomeSequenceId];
          bestStrand    = states->sw->strands[genomeSequenceId];
          memcpy(bestStates, states->sw->states, alignmentSize * sizeof(state_t));
        }
      }
    }
  }
  if (bestNErrors > parameters->maxNErrors) {
    path->nCells  = previousNCells;
    path->depth   = previousDepth;
    path->readPos = previousReadPos;
    //printf("    Leaving tryShortCuts wrong with cell %" PRIu64 " at depth %zu, last nt %i, and read pos %zu\n", path->cellIds[path->nCells], path->depth, path->nucleotides[path->depth-1], path->readPos);
    return false;
  }
  assert(alignmentSize > 0);
  previousState = &states->states[previousDepth][baseNErrors][bestStateId];
  //printf("\tSize of alignment: %u, starting from %p: %lu, %u, %u\n", alignmentSize, previousState, previousDepth, baseNErrors, bestStateId);
  currentNErrors = baseNErrors;
  genomeAlignmentSize = 0;
  if (bestNErrors == 0) {
    outputSam->backtraceSize       = 1;
    outputSam->backtraceLengths[0] = path->depth;
    outputSam->backtraceCigar[0]   = 'M';
  }
  else {
    computeBacktrace(previousState, outputSam);
  }
  for (unsigned int i = 0; i < alignmentSize; ++i) {
    nextState = addState(states, path->depth, currentNErrors);
    if (nextState == NULL) {
      return false;
    }
    setState(nextState, &bestStates[i].interval, bestStates[i].trace, previousState);
    if ((nextState->trace & BACKTRACE_MASK) != DELETION) {
      //printf("\t\tAdding nucleotide %c\n", "ACGT"[bestStates[i].trace & NUCLEOTIDE_MASK]);
      appendNucleotidePath(path, bestStates[i].trace & NUCLEOTIDE_MASK, DNA5_TO_CHAR[bestStates[i].trace & NUCLEOTIDE_MASK]);
    }
    if (! hasTrace(nextState, MATCH)) {
      ++currentNErrors;
    }
    if (! hasTrace(nextState, INSERTION)) {
      ++genomeAlignmentSize;
    }
    //printf("\tAdding state @ %lu/%u/%lu: ", path->depth, currentNErrors, states->nStates[path->depth][currentNErrors]-1);
    //printState(&states->states[path->depth][currentNErrors][states->nStates[path->depth][currentNErrors]-1], 100);
    previousState = nextState;
    assert(currentNErrors <= bestNErrors);
    cigar = CIGAR[nextState->trace >> BACKTRACE_OFFSET];
    if (bestNErrors != 0) {
      if ((outputSam->backtraceSize == 0) || (outputSam->backtraceCigar[outputSam->backtraceSize-1] != cigar)) {
        outputSam->backtraceCigar[outputSam->backtraceSize] = cigar;
        outputSam->backtraceLengths[outputSam->backtraceSize] = 1;
        ++outputSam->backtraceSize;
      }
      else {
        ++outputSam->backtraceLengths[outputSam->backtraceSize-1];
      }
    }
  }
  assert(currentNErrors == bestNErrors);
  //printStates(states, path->depth+1);
  //printPath(path);
  writeQname(outputSam, cell->counts);
  memcpy(outputSam->forwardSeq,  path->read + path->readPos,  (path->depth+1) * sizeof(char));
  memcpy(outputSam->forwardQual, findQuality(&tree->qualities, cellId), (path->depth+1) * sizeof(char));
  outputSam->isBackwardSet = false;
  computeCigar(outputSam);
  bestPos = (bestStrand)? bestPos + genomeAlignmentSize: bestPos - genomeAlignmentSize;
  rid     = bns_pos2rid(bns, bestPos);
  //printf("Seq: %s, qual: %s, strand: %i, cigar F: %s, cigar B: %s, qname: %s\n", outputSam->forwardSeq, outputSam->forwardQual, bestStrand, outputSam->forwardCigar, outputSam->backwardCigar, outputSam->qname);
  preparePrintRead(bestPos, rid, bestStrand, 1, 1, bestNErrors, outputSam);
  //printf("    Leaving tryShortCuts right with cell %" PRIu64 " at depth %zu, last nt %i, and read pos %zu\n", path->cellIds[path->nCells], path->depth, path->nucleotides[path->depth-1], path->readPos);
  //printPath(path);
  path->nCells  = previousNCells;
  path->depth   = previousDepth;
  path->readPos = previousReadPos;
  backtrackStates(states, path->depth);
  return true;
}

/*
bool tryShortCut (const tree_t *tree, states_t *states, path_t *path, FILE *outputSamFile) {
  return false;
  shortCut_t *shortCut = path->shortCut;
  size_t depth = path->depth;
  uint64_t cellId = path->cellIds[depth];
  unsigned short nextNucleotide;
  cell_t *cell = &tree->cells[cellId];
  bool lastNucleotide;
  size_t readPos = path->readPos;
  char *quality;
  //printf("Trying short cut.  Min errors: %zu, # states: %zu, depth: %zu\n", states->minErrors[depth], states->nStates[depth][0], depth);
  //printf("  k: %" PRIu64", l: %" PRIu64 "\n", states->states[path->depth][0][0].k, states->states[path->depth][0][0].l);
  //printPath(path);
  if (! shortCut->isSet) {
    setShortCut(shortCut, getStateInterval(&states->states[path->depth][0][0])->k, depth);
  }
  else {
    resetShortCut(shortCut, depth);
  }
  //printf("  Step 1 ok, base nt is '%c'\n", DNA5_TO_CHAR_REV[is_rev][_get_pac(pac, pos)]);
  ++stats->nShortCuts;
  while (true) {
    //printf("  Depth %zu\n", path->depth);
    incShortCut(shortCut);
    nextNucleotide = DNA5_TO_INT_REV[shortCut->isRev][_get_pac(pac, shortCut->pos)];
    //printf("  Next nucleotide: %c\n", DNA5_TO_CHAR[nextNucleotide]);
    cellId = cell->children[nextNucleotide];
    lastNucleotide = (cellId == NO_DATA);
    for (unsigned short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
      if ((nucleotide != nextNucleotide) && (cell->children[nucleotide] != NO_DATA)) {
        //printf("    Other nucleotide: %i is used, exiting at depth %zu/%zu.\n", nucleotide, path->depth, shortCut->depth);
        addMiss(shortCut);
        return true;
      }
    }
    if (lastNucleotide) {
      //printf("    Last nucleotide @ depth %zu/%zu, exiting.\n", path->depth, shortCut->depth);
      ++stats->nShortCutSuccesses;
      return (goRightTree(tree, path));
    }
    --readPos;
    path->read[readPos] = DNA5_TO_CHAR[nextNucleotide];
    cell = &tree->cells[cellId];
    if ((quality = findQuality(&tree->qualities, cellId)) != NULL) {
      //printf("    Printing read @ depth %zu\n", path->depth);
      printReadUniqueNoError(shortCut->isRev? 0: 1, bns->anns[shortCut->rid].name, shortCut->pos, shortCut->depth, path->read + readPos, quality, cell->counts, outputSamFile);
    }
  }
  assert(false);
  return true;
}
*/

void printProgress(path_t *path) {
  if (path->depth == TREE_BASE_SIZE) {
    if ((path->cellIds[path->nCells] & 11111) == 0) {
      fprintf(stderr, "Progress: %" PRIu64 "/%u\n", path->cellIds[path->nCells], N_TREE_BASE);
    }
  }
}

/**
 * Map one nucleotide
 * Here, we have added a new nucleotide in (reads) path.
 * We are looking for the corresponding mappings at the path->depth level.
 */
bool findBestMapping (states_t *states, path_t *path) {
  //printf("Finding best mapping\n");
  //printPath(path);
  //printStates(states, path->depth);
  //printf("\t\t\t\tRead: %s\n", path->read + path->readPos);
  //for (size_t i = 0; i < states->nStates[path->depth-1][states->minErrors[path->depth-1]]; ++i) { printState(&states->states[path->depth-1][states->minErrors[path->depth-1]][i], path->maxDepth); }
  printProgress(path);
  /*
  if (path->depth <= TREE_BASE_SIZE) {
    //printf("  with errors\n");
    mapWithErrors(states, path);
    return true;
  }
  */
  if (mapWithoutError(states, path->depth, path->nucleotides[path->depth-1], states->minErrors[path->depth-1])) {
    //printf("    ... without error\n");
    return true;
  }
  //printf("    # errors: %zu/%zu\n", states->minErrors[path->depth-1], parameters->maxNErrors);
  if (states->minErrors[path->depth-1] >= parameters->maxNErrors) {
    //printf("    ... exiting\n");
    return false;
  }
  //printf("    ... adding error\n");
  return addError(states, path);
}

/**
 * Do the mapping
 */
void _map (const tree_t *tree, states_t *states, path_t *path, outputSam_t *outputSam) {
  bool mappable = true;
  char *quality;
  while (true) {
    if ((mappable) && (shortCutCondition(states, tree, path))) {
      tryShortCuts(tree, states, path, outputSam);
      mappable = false;
      /*
      if (mappable) {
        //printf("Short cut with positive exit\n"); fflush(stdout);
        if ((quality = findQuality(&tree->qualities, path->cellIds[path->nCells])) != NULL) {
          printRead(states, path, quality, tree->cells[path->cellIds[path->nCells]].counts, outputSam);
        }
      }
      */
    }
    else {
      // advance in the tree
      //printf("Going to next tree\n");
      if (! goNextTree(tree, states, path, mappable)) {
        return;
      }
      //printPath(path); fflush(stdout);
      //printf("%" PRIu64 ": ", path->cellIds[path->nCells]);
      //printCell(&tree->cells[path->cellIds[path->nCells]]);
      //printf("\n");
      // prefix of the read is not mappable
      mappable = findBestMapping(states, path);
      //printf("\tRead is mappable: %s\n", (mappable)? "yes": "no");
      if ((path->depth >= TREE_BASE_SIZE) && (mappable) && (path->edgeLength == 0)) {
        if ((quality = findQuality(&tree->qualities, path->cellIds[path->nCells])) != NULL) {
          printRead(states, path, quality, tree->cells[path->cellIds[path->nCells]].counts, outputSam);
        }
      }
    }
  }
}

/**
 * Allocate/free structures before/after mapping
 */
void map (const tree_t *tree, outputSam_t *outputSam) {
  states_t *states = initializeStates(tree->depth);
  path_t   *path   = initializePath(tree->depth);
  //printf("depth: %zu\n", tree->depth);
  //addState(&states, 0, 0, firstState);
  //goDownTree(tree, &path);
  _map(tree, states, path, outputSam);
  freeStates(states);
  freePath(path);
}

FILE *openSamFile() {
  FILE *outputSamFile = fopen(parameters->outputSamFileName, "w");
  if (outputSamFile == NULL) {
    return NULL;
  }
  fprintf(outputSamFile, "@HD\tVN:1.6\tSO:unsorted\n");
  for (int i = 0; i < bns->n_seqs; ++i) {
    fprintf(outputSamFile, "@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
  }
  return outputSamFile;
}

int main(int argc, char const ** argv) {
  int returnCode = 0;
  parameters_t param;
  stats_t stat;
  tree_t tree;
  outputSam_t outputSam;
  bwaidx_t *idx = NULL;
  nReads = 0;
  parameters = &param;
  stats = &stat;
  initializeStats();
  returnCode = parseCommandLine(argc, argv);
  if (returnCode != EXIT_SUCCESS) return returnCode;
  createTree(&tree);
  for (unsigned int fileId = 0; fileId < parameters->nReadsFiles; ++fileId) {
    printf("Reading reads file %s...\n", parameters->readsFileNames[fileId]);
    returnCode = readReadsFile(parameters->readsFileNames[fileId], &tree, fileId);
    if (returnCode != EXIT_SUCCESS) return returnCode;
    puts("... done.");
  }
  puts("Filtering tree...");
  filterTree(&tree);
  printf("... done.\n");
  printf("Maximum read size: %zu\n", tree.depth);
  computeTreeStats(&tree);
  if (parameters->outputReadsFileName != NULL) {
    puts("Printing tree...");
    returnCode = printTree(parameters->outputReadsFileName, &tree);
    if (returnCode != EXIT_SUCCESS) return returnCode;
    puts("... done.");
  }
  puts("Loading genome...");
  idx = loadGenomeFile(parameters->genomeFileName);
  if (idx == NULL) return EXIT_FAILURE;
  puts("... done.");
  puts("Mapping...");
  pac = idx->pac;
  bwt = idx->bwt;
  bns = idx->bns;
  FILE *outputSamFile = openSamFile();
  if (outputSamFile == NULL) {
    printf("Error!  Cannot write to output SAM file '%s'.\nExiting.\n", param.outputSamFileName);
  }
  outputSam.file = outputSamFile;
  createOutputSam(&outputSam, tree.depth);
  map(&tree, &outputSam);
  freeTree(&tree);
  freeOutputSam(&outputSam);
  bwa_idx_destroy(idx);
  puts("... done.");
  printStats();
  return 0;
}
