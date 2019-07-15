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

#define TREE_BASE_SIZE   8
#define N_TREE_BASE  65536

#define N_STATES        0x5000

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
#define BACKTRACE_MASK    3
#define PREPROCESSED     16

#define CIGAR_SECONDARY_HIT 0x100
#define CIGAR_REVERSE 0x10

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

char *reverseSequence (char *sequence, size_t size) {
  char *revSequence = (char *) malloc(sizeof(char) * (size+1));
  for (size_t i = 0; i < size; ++i) {
    revSequence[i] = sequence[size-i-1];
  }
  revSequence[size] = 0;
  return revSequence;
}

char *reverseComplementSequence (char *sequence, size_t size) {
  char *revSequence = (char *) malloc(sizeof(char) * (size+1));
  for (size_t i = 0; i < size; ++i) {
    revSequence[i] = COMP[(int) sequence[size-i-1]];
  }
  revSequence[size] = 0;
  return revSequence;
}

typedef struct {
  unsigned long int nReads;
  unsigned long int nShortReads;
  unsigned long int nShortCuts;
  unsigned long int nShortCutSuccesses;
  unsigned long int nDown;
  unsigned long int nDownPreprocessed;
  size_t            maxNStates;
} stats_t;

stats_t *stats;

void initializeStats () {
  stats->nReads             = 0;
  stats->nShortReads        = 0;
  stats->nShortCuts         = 0;
  stats->nShortCutSuccesses = 0;
  stats->nDown              = 0;
  stats->nDownPreprocessed  = 0;
  stats->maxNStates         = 0;
}

void printStats () {
  char *savedLocale;
  savedLocale = setlocale (LC_ALL, NULL);
  setlocale(LC_NUMERIC, "");
  printf("Very small sequences: %'lu/%'lu\n", stats->nShortReads, stats->nReads);
  printf("Short cut succeeded %'lu/%'lu\n", stats->nShortCutSuccesses, stats->nShortCuts);
  printf("# BWT preprocessed %'lu/%'lu\n", stats->nDownPreprocessed, stats->nDown);
  printf("# max states %'zu/%'i\n", stats->maxNStates, N_STATES);
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
  uint64_t children [N_NUCLEOTIDES];
  count_t *counts;
} cell_t;

void createCell (cell_t *cell) {
  for (unsigned short i = 0; i < N_NUCLEOTIDES; ++i) {
    cell->children[i] = NO_DATA;
  }
  //printf("%u %zu %p\n", parameters->nReadsFiles, sizeof(count_t), cell);
  cell->counts = (count_t *) calloc(parameters->nReadsFiles, sizeof(count_t));
}

void freeCell (cell_t *cell) {
  free(cell->counts);
}

typedef struct {
  char       **qualities;
  //uint64_t    *cellIds;
  unsigned int nQualities;
  unsigned int nAllocated;
} quality_t;

void createQualities (quality_t *qualities) {
  qualities->nAllocated = INIT_N_QUALITIES;
  qualities->qualities  = (char **) calloc(qualities->nAllocated, sizeof(char *));
  qualities->nQualities = 0;
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
  if ((cellId >= qualities->nQualities) || (qualities->qualities[cellId] == NULL)) return NO_QUALITY;
  return cellId;
}

char *findQuality (const quality_t *qualities, uint64_t cellId) {
  unsigned int qualityId = findQualityId(qualities, cellId);
  if (qualityId == NO_QUALITY) return NULL;
  return qualities->qualities[qualityId]; 
}

void addQuality (quality_t *qualities, uint64_t cellId, size_t l, char *quality) {
  while (cellId >= qualities->nAllocated) {
    qualities->nAllocated *= 2;
  }
  if ((qualities->qualities = (char **) realloc(qualities->qualities, qualities->nAllocated * sizeof(char *))) == NULL) {
    printf("Cannot allocate memory for qualities of size %u.\nExiting.\n", qualities->nAllocated);
    exit(EXIT_FAILURE);
  }
  qualities->qualities[cellId] = strndup(quality, l);
}

void replaceQuality (quality_t *qualities, unsigned int qualityId, size_t l, char *quality) {
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
  unsigned int c = 0;
  for (short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    if (cell->children[nucleotide] != NO_DATA) {
      ++c;
      ++(*nNodes);
      if (findQualityId(&tree->qualities, cell->children[nucleotide]) != NO_QUALITY) {
        ++(*nQualities);
      }
    }
  }
  if (c == 1) {
    ++branchSize;
  }
  else {
    ++branchSizes[branchSize];
    branchSize = 0;
  }
  for (short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    if (cell->children[nucleotide] != NO_DATA) {
      _computeTreeStats(tree, stats, statsSum, branchSizes, branchSize, &tree->cells[cell->children[nucleotide]], depth+1, nNodes, nQualities);
    }
  }
  ++stats[depth][c];
  ++statsSum[c];
}

void computeTreeStats (const tree_t *tree) {
  unsigned int **stats = (unsigned int **) malloc((tree->depth+1) * N_NUCLEOTIDES * sizeof(unsigned int *));
  unsigned int statsSum[N_NUCLEOTIDES] = { 0, 0, 0, 0 };
  unsigned int *branchSizes = (unsigned int *) calloc(tree->depth+1, sizeof(unsigned int));
  unsigned int nNodes = 0;
  unsigned int nQualities = 0;
  unsigned int s;
  for (size_t depth = 0; depth <= tree->depth; ++depth) {
    stats[depth] = (unsigned int *) calloc(N_NUCLEOTIDES, sizeof(unsigned int));
  }
  for (size_t cellId = 0; cellId < N_TREE_BASE; ++cellId) {
    _computeTreeStats(tree, stats, statsSum, branchSizes, 0, &tree->cells[cellId], TREE_BASE_SIZE-1, &nNodes, &nQualities);
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

uint64_t goDown (tree_t *tree, uint64_t cellId, unsigned short childId) {
  cell_t *cell = &tree->cells[cellId];
  uint64_t newCellId = cell->children[childId];
  if (newCellId != NO_DATA) {
    return newCellId;
  }
  newCellId = addCell(tree);
  cell = &tree->cells[cellId]; // may be reallocated!
  cell->children[childId] = newCellId;
  return newCellId;
}

bool addSequence (tree_t *tree, size_t l, char *sequence, char *quality, unsigned int fileId) {
  uint64_t cellId = 0;
  int sequenceId;
  assert(strlen(sequence) == strlen(quality));
  assert(strlen(quality) == l);
  for (int i = 0; i < TREE_BASE_SIZE; ++i) {
    sequenceId = l - i - 1;
    if (sequenceId < 0) {
      return false;
    }
    cellId <<= NUCLEOTIDES_BITS;
    cellId += CHAR_TO_DNA5[(int) sequence[sequenceId]];
  }
  printf("First id: %lu, %s\n", cellId, sequence);
  for (--sequenceId; sequenceId >= 0; --sequenceId) {
    cellId = goDown(tree, cellId, CHAR_TO_DNA5[(int) sequence[sequenceId]]);
  }
  setQuality(tree, cellId, l, quality, fileId);
  tree->depth = MAX(tree->depth, l);
  return true;
}

void __printTree (const tree_t *tree, FILE *outFile, uint64_t *readId, char *read, size_t readPos, uint64_t cellId) {
  uint64_t nextCellId;
  cell_t *cell = &tree->cells[cellId];
  char *quality;
  //read[readPos] = 0; printf("Base triplet: %d, size: %zu, read: %s\n", triplet, readPos, read);
  //read[readPos] = 0; printf("Current state: %zu, %s, %zu, %lu\n", *readId, read, readPos, cellId); fflush(stdout);
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
  //printf("\tNew triplet: %d\n", triplet);
  for (unsigned short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    //printf("\tTrying child %d\n", nucleotide); fflush(stdout);
    nextCellId = cell->children[nucleotide];
    if (nextCellId != NO_DATA) {
      //printf("\tNext triplet (%d): %d\n", i, nextTriplet);
      read[tree->depth-readPos-1] = DNA5_TO_CHAR[nucleotide];
      __printTree(tree, outFile, readId, read, readPos+1, nextCellId);
    }
  }
}

void _printTree (const tree_t *tree, FILE *outFile, uint64_t *readId, char *read, size_t readPos, uint64_t cellId) {
  if (readPos == TREE_BASE_SIZE) {
    __printTree (tree, outFile, readId, read, readPos, cellId);
    return;
  }
  cellId <<= NUCLEOTIDES_BITS;
  for (unsigned short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    read[tree->depth-readPos-1] = DNA5_TO_CHAR[nucleotide];
    _printTree(tree, outFile, readId, read, readPos+1, cellId+nucleotide);
  }
}

int printTree (char *fileName, const tree_t *tree) {
  FILE *outFile = fopen(fileName, "w");
  if (outFile == NULL) return EXIT_FAILURE;
  uint64_t readId = 0;
  char *read = (char *) malloc((tree->depth+1) * sizeof(char));
  read[tree->depth] = 0;
  _printTree(tree, outFile, &readId, read, 0, 0);
  free(read);
  fclose(outFile);
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

bool __filterTree (const tree_t *tree, size_t readPos, uint64_t cellId, unsigned short triplet, count_t *tripletCount) {
  uint64_t nextCellId;
  unsigned short nextTriplet;
  bool foundRead = false;
  cell_t *cell = &tree->cells[cellId];
  triplet &= TRIPLET_MASK;
  triplet <<= NUCLEOTIDES_BITS;
  if (findQuality(&tree->qualities, cellId) != NULL) {
    foundRead = true;
  }
  for (unsigned short i = 0; i < N_NUCLEOTIDES; ++i) {
    nextCellId = cell->children[i];
    if (nextCellId != NO_DATA) {
      nextTriplet = triplet | i;
      if (readPos >= TRIPLET-1) tripletCount[nextTriplet] += 1;
      if (tripletCount[nextTriplet] > parameters->lowComplexityThreshold) {
        cell->children[i] = NO_DATA;
      }
      else {
        if (! __filterTree(tree, readPos+1, nextCellId, nextTriplet, tripletCount)) {
          cell->children[i] = NO_DATA;
        }
        else {
          foundRead = true;
        }
      }
      if (readPos >= TRIPLET-1) tripletCount[nextTriplet] -= 1;
    }
  }
  return foundRead;
}

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
    if (readPos >= TRIPLET-1) tripletCount[nextTriplet] += 1;
    if (__filterTree(tree, readPos+1, cellId + nucleotide, nextTriplet, tripletCount)) {
      foundRead = true;
    }
    if (readPos >= TRIPLET-1) tripletCount[nextTriplet] -= 1;
  }
  return foundRead;
}

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

typedef struct {
  bwtint_t k, l;
} bwtinterval_t;

typedef struct preprocessedInterval_t preprocessedInterval_t;

struct preprocessedInterval_t {
  bwtinterval_t interval;
};

typedef struct {
  preprocessedInterval_t *intervals;
  size_t nIntervals;
} preprocessedIntervals_t;

preprocessedIntervals_t *initializePreprocessedStates () {
  preprocessedIntervals_t *preprocessedIntervals = (preprocessedIntervals_t *) malloc(sizeof(preprocessedIntervals_t));
  preprocessedIntervals->nIntervals = (powl(N_NUCLEOTIDES, PREPROCESSED_DEPTH+1) - 1) / (N_NUCLEOTIDES - 1);
  //printf("#intervals: %zu\n", preprocessedIntervals->nIntervals);
  preprocessedIntervals->intervals = (preprocessedInterval_t *) malloc(preprocessedIntervals->nIntervals * sizeof(preprocessedInterval_t));
  return preprocessedIntervals;
}

void freePreprocessedStates (preprocessedIntervals_t *preprocessedIntervals) {
  free(preprocessedIntervals->intervals);
}

//static unsigned long nIntervals;

void _populatePreprocessedIntervals (preprocessedIntervals_t *preprocessedIntervals, size_t depth, size_t id) {
  size_t nextId = id * N_NUCLEOTIDES + 1;
  ++depth;
  for (unsigned int nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    bwt_2occ(bwt, preprocessedIntervals->intervals[id].interval.k, preprocessedIntervals->intervals[id].interval.l, nucleotide, &preprocessedIntervals->intervals[nextId].interval.k, &preprocessedIntervals->intervals[nextId].interval.l);
    preprocessedIntervals->intervals[nextId].interval.k = bwt->L2[nucleotide] + preprocessedIntervals->intervals[nextId].interval.k + 1;
    preprocessedIntervals->intervals[nextId].interval.l = bwt->L2[nucleotide] + preprocessedIntervals->intervals[nextId].interval.l;
    if ((preprocessedIntervals->intervals[nextId].interval.k <= preprocessedIntervals->intervals[nextId].interval.l) && (depth < PREPROCESSED_DEPTH)) {
      _populatePreprocessedIntervals(preprocessedIntervals, depth, nextId);
      //printf("id: %zu\n", nextId);
    }
    ++nextId;
    //printf("  %c: %" PRIu64 "-%" PRIu64 "...\n", "ACGT"[nucleotide], preprocessedFreeState->interval.k, preprocessedFreeState->interval.l);
  }
}

void populatePreprocessedIntervals (preprocessedIntervals_t *preprocessedIntervals) {
  //nIntervals = 0;
  printf("Pre-processing %zu intervals...\n", preprocessedIntervals->nIntervals);
  preprocessedIntervals->intervals[0].interval.k = 0;
  preprocessedIntervals->intervals[0].interval.l = bwt->bwt_size;
  _populatePreprocessedIntervals(preprocessedIntervals, 0, 0);
  //printf("#intervals used: %lu\n", nIntervals);
  puts("... done.");
}

typedef struct state_t state_t;

struct state_t {
  union {
    bwtinterval_t interval;
    preprocessedInterval_t *preprocessedInterval;
  };
  unsigned char trace;
  struct state_t *previousState;
};

bwtinterval_t *getStateInterval (state_t *state) {
  if (state->trace & PREPROCESSED) return &state->preprocessedInterval->interval;
  return &state->interval;
}

void printState(state_t *state, size_t maxDepth) {
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
  printf("\t\t\t\t%" PRIu64 "-%" PRIu64 ": %s (%c/%c -> %p)\n", getStateInterval(state)->k, getStateInterval(state)->l, tmpSeq2, CIGAR[(state->trace >> BACKTRACE_OFFSET) & BACKTRACE_MASK], "ACGT"[state->trace & NUCLEOTIDE_MASK], state->previousState);
}

typedef struct {
  state_t ***states;
  size_t **nStates;
  size_t *nStatesPerPosition;
  size_t *minErrors;
  size_t *maxErrors;
  size_t treeSize;
  preprocessedIntervals_t *preprocessedIntervals;
} states_t;

bool goDownBwt (const states_t *states, state_t *previousState, unsigned short nucleotide, state_t *newState) {
  ++stats->nDown;
  //printf("    Going down BWT from range %" PRId64 "-%" PRIu64 " and nt %hu, preprocessed: %s\n", getStateInterval(previousState)->k, getStateInterval(previousState)->l, nucleotide, (previousState->trace & PREPROCESSED)? "true": "false");
  if (previousState->trace & PREPROCESSED) {
    size_t nextId = (previousState->preprocessedInterval - states->preprocessedIntervals->intervals) * N_NUCLEOTIDES + nucleotide + 1;
    //printf("      preprocessed next id should be %zu -> %zu (%" PRIu64 "-%" PRIu64 ")\n", previousState->preprocessedInterval - states->preprocessedIntervals->intervals, nextId, states->preprocessedIntervals->intervals[nextId].interval.k, states->preprocessedIntervals->intervals[nextId].interval.l);
    if (nextId < states->preprocessedIntervals->nIntervals) {
      //printf("      ok\n");
      ++stats->nDownPreprocessed;
      //printf("      first case\n");
      newState->preprocessedInterval = &states->preprocessedIntervals->intervals[nextId];
      newState->trace                = PREPROCESSED;
      //if (newState->preprocessedInterval->interval.k <= newState->preprocessedInterval->interval.l) printf("      ok!\n");
      return (newState->preprocessedInterval->interval.k <= newState->preprocessedInterval->interval.l);
    }
  }
  bwt_2occ(bwt, getStateInterval(previousState)->k-1, getStateInterval(previousState)->l, nucleotide, &newState->interval.k, &newState->interval.l);
  newState->interval.k = bwt->L2[nucleotide] + newState->interval.k + 1;
  newState->interval.l = bwt->L2[nucleotide] + newState->interval.l;
  //if (newState->interval.k <= newState->interval.l) printf("      ok!\n");
  newState->trace = 0;
  return (newState->interval.k <= newState->interval.l);
}

void printStates (states_t *states, size_t depth) {
  for (size_t i = 0; i <= depth; ++i) {
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

bool compareStates (state_t *state1, state_t *state2) {
  return ((getStateInterval(state1)->k == getStateInterval(state2)->k) && (getStateInterval(state1)->l == getStateInterval(state2)->l));
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

size_t simplifyStates (state_t *states, size_t nStates) {
  //printf("\t\t\t\tSimplify states from %zu ", nStates);
  qsort(states, nStates, sizeof(state_t), sortCompareStates);
  size_t firstStateId = 0;
  for (size_t secondStateId = 1; secondStateId < nStates; ++secondStateId) {
    if (! compareStates(&states[firstStateId], &states[secondStateId])) {
    //if (! canMerge(&states[firstStateId], &states[secondStateId])) {
      ++firstStateId;
      states[firstStateId] = states[secondStateId];
    }
  }
  //printf("to %zu\n", firstStateId+1);
  return (firstStateId+1);
}

bool addState(states_t *states, size_t depth, size_t nErrors, state_t *state) {
  //printf("\t\t\tAdding one state (%" PRIu64 ", %" PRIu64 ") at (depth = %zu, # errors = %zu), %zu/%d occupied\n", state->k, state->l, depth, nErrors, states->nStates[depth][nErrors], N_STATES);
  //printState(state, states->treeSize);
  if (states->nStates[depth][nErrors] == N_STATES-1) {
    //printf("Exiting because # states = %zu >= %i at depth %zu with %zu errors, before simplification.\n", states->nStates[depth][nErrors], N_STATES, depth, nErrors);
    //printStates(states, depth);
    states->nStates[depth][nErrors] = N_STATES;
    return false;
  }
  states->states[depth][nErrors][states->nStates[depth][nErrors]] = *state;
  ++states->nStates[depth][nErrors];
  ++states->nStatesPerPosition[depth];
  stats->maxNStates = MAX(stats->maxNStates, states->nStates[depth][nErrors]);
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
  return true;
}

states_t *initializeStates(size_t treeSize) {
  states_t *states           = (states_t *) malloc(sizeof(states_t));
  states->treeSize           = treeSize;
  states->states             = (state_t ***) malloc((treeSize+1) * sizeof(states_t **));
  states->nStates            = (size_t **)   malloc((treeSize+1) * sizeof(size_t *));
  states->nStatesPerPosition = (size_t *)    calloc((treeSize+1),  sizeof(size_t));
  states->minErrors          = (size_t *)    malloc((treeSize+1) * sizeof(size_t));
  states->maxErrors          = (size_t *)    malloc((treeSize+1) * sizeof(size_t));
  for (size_t depth = 0; depth <= treeSize; ++depth) {
    states->states[depth]    = (state_t **) malloc((parameters->maxNErrors+1) * sizeof(states_t *));
    states->nStates[depth]   = (size_t *)   calloc((parameters->maxNErrors+1),  sizeof(size_t));
    states->minErrors[depth] = SIZE_MAX;
    states->maxErrors[depth] = SIZE_MAX;
    for (size_t nErrors = 0; nErrors <= parameters->maxNErrors; ++nErrors) {
      states->states[depth][nErrors] = (state_t *) malloc(N_STATES * sizeof(states_t));
    }
  }
  preprocessedIntervals_t *preprocessedIntervals = initializePreprocessedStates();
  populatePreprocessedIntervals(preprocessedIntervals);
  states->preprocessedIntervals = preprocessedIntervals;
  state_t baseState;
  baseState.preprocessedInterval = &preprocessedIntervals->intervals[0];
  baseState.trace                = PREPROCESSED;
  baseState.previousState        = NULL;
  addState(states, 0, 0, &baseState);
  return states;
}

void backtrackStates(states_t *states, size_t level) {
  //printf("\t\t\tBacktracking to level %zu\n", level);
  for (size_t i = level; i < states->treeSize; ++i) {
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
  for (size_t i = 0; i < states->treeSize; ++i) {
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
  freePreprocessedStates(states->preprocessedIntervals);
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

typedef struct {
  unsigned short *nucleotides;
  uint64_t       *cellIds;
  size_t          maxDepth;
  size_t          depth;
  char           *read;
  size_t          readPos;
  shortCut_t     *shortCut;
} path_t;

void printPath (path_t *path) {
  assert(path->depth < path->maxDepth);
  printf("\t\t\t\t");
  for (size_t i = 0; i < path->depth; ++i) {
    assert(path->nucleotides[i] < N_NUCLEOTIDES);
    printf("%c ", "ACGT"[(int) path->nucleotides[i]]);
  }
  printf("\n");
}

path_t *initializePath (size_t maxDepth) {
  path_t *path         = (path_t *)         malloc(sizeof(path_t));
  path->nucleotides    = (unsigned short *) malloc(maxDepth * sizeof(unsigned short));
  path->cellIds        = (uint64_t *)       malloc((maxDepth+1) * sizeof(uint64_t));
  path->read           = (char *)           malloc((maxDepth+1) * sizeof(char));
  path->maxDepth       = maxDepth;
  path->depth          = 0;
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

void appendPath (path_t *path, size_t depth, uint64_t cellId, short nucleotide, char c) {
  path->nucleotides[depth] = nucleotide;
  ++path->depth;
  --path->readPos;
  path->cellIds[path->depth] = cellId;
  path->read[path->readPos] = c;
}

bool goDownTree (const tree_t *tree, path_t *path) {
  uint64_t cellId;
  //printf("  Go down from height %zu with read '%s' ", path->depth, path->read+path->readPos);
  if (path->depth < TREE_BASE_SIZE) {
    cellId = path->cellIds[path->depth];
    cellId <<= NUCLEOTIDES_BITS;
    appendPath(path, path->depth, cellId, 0, 'A');
    return true;
  }
  cell_t *cell = &tree->cells[path->cellIds[path->depth]];
  for (unsigned short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    if ((cellId = cell->children[nucleotide]) != NO_DATA) {
      appendPath(path, path->depth, cellId, nucleotide, DNA5_TO_CHAR[nucleotide]);
      //printf("to %zu with read '%s'\n", path->depth, path->read+path->readPos);
      return true;
    }
  }
  //printf("to nothing\n");
  return false;
}

bool goRightTreeBase (const tree_t *tree, path_t *path) {
  uint64_t cellId, nextCellId;
  assert(path->depth < TREE_BASE_SIZE);
  if (path->depth == 0) {
    return false;
  }
  cellId = path->nucleotides[tree->depth];
  nextCellId = cellId + 1;
  while (true) {
    if (cellId == nextCellId) {
      return true;
    }
    path->cellIds[path->depth] = nextCellId;
    path->read[path->readPos] = DNA5_TO_CHAR[nextCellId & NUCLEOTIDE_MASK];
    path->nucleotides[path->depth] = nextCellId & NUCLEOTIDE_MASK; 
    if (path->depth == 0) {
      return false;
    }
    --path->depth;
    ++path->readPos;
    cellId     >>= NUCLEOTIDES_BITS;
    nextCellId >>= NUCLEOTIDES_BITS;
  }
  return false;
}

bool goRightTreeNotBase (const tree_t *tree, path_t *path) {
  cell_t *cell;
  //printf("  Go right read from depth %zu with read '%s'\n", path->depth, path->read+path->readPos);
  while (true) {
    --path->depth;
    ++path->readPos;
    if (path->depth < TREE_BASE_SIZE) {
      return goRightTreeBase(tree, path);
    }
    //printf("    ... trying depth %zu\n", path->depth);
    for (size_t nucleotide = path->nucleotides[path->depth]+1; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
      cell = &tree->cells[path->cellIds[path->depth]];
      //printf("      ... nucleotide is %zu, cell is %" PRIu64 "\n", nucleotide, tree->cells[path->cellIds[path->depth]].children[nucleotide]); 
      if ((cellId = cell->children[nucleotide]) != NO_DATA) {
        path->nucleotides[path->depth] = nucleotide;
        ++path->depth;
        --path->readPos;
        path->cellIds[path->depth] = cellId;
        path->read[path->readPos] = DNA5_TO_CHAR[nucleotide];
        //printf("    ... going to %zu with read '%s'\n", path->depth, path->read+path->readPos);
        return true;
      }
    }
  }
  //printf("    ... going to nothing\n");
  return false;
}

bool goRightTree (const tree_t *tree, path_t *path) {
  unsetShortCut(path->shortCut);
  if (path->depth < TREE_BASE_SIZE) {
    return goRightTreeNotBase(tree, path);
  }
  return goRightTreeBase(tree, path);
}

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
  if (! goRightTree(tree, path)) {
    return false;
  }
  //printf("  Now\n");
  //printStates(states, path->depth);
  backtrackStates(states, path->depth);
  return true;
}

void writeQname (char *qname, count_t *counts) {
  size_t qnameLength = sprintf(qname, "read%lu_x", ++nReads);
  for (size_t readsFileId = 0; readsFileId < parameters->nReadsFiles; ++readsFileId) {
    qnameLength += sprintf(qname+qnameLength, "_%lu", counts[readsFileId]);
  }
}

void printReadUniqueNoError (int strand, char *chrName, int64_t pos, size_t readLength, char *sequence, char *quality, count_t *counts, FILE *outputSamFile) {
  char qname[255];
  unsigned int flag = 0;
  char *seq = sequence;
  char *qual = quality;
  if (strand == 0) {
    flag = CIGAR_REVERSE;
    seq  = reverseComplementSequence(sequence, readLength);
    qual = reverseSequence(quality, readLength);
  }
  writeQname(qname, counts);
  fprintf(outputSamFile, "%s\t%u\t%s\t%" PRId64 "\t40\t%zuM\t*\t0\t0\t%s\t%s\tNH:i:1\tHI:i:1\tIH:i:1\tNM:i:0\n", qname, flag, chrName, pos, readLength, seq, qual);
  if (strand == 0) {
    free(seq);
    free(qual);
  }
}

void printRead (states_t *states, path_t *path, char *quality, count_t *counts, FILE *outputSamFile) {
  size_t depth = path->depth;
  size_t readLength = depth;
  char qname[255];
  unsigned int flag;
  int64_t pos;
  unsigned int mapq;
  char *cigar;
  unsigned char backtraceLengths[255];
  char backtraceCigar[255];
  size_t backtraceSize;
  char forwardCigar[255];
  char backwardCigar[255];
  size_t forwardCigarSize = 0;
  size_t backwardCigarSize = 0;
  char *forwardSeq = path->read + path->readPos;
  char *backwardSeq = NULL;
  char *seq;
  char *forwardQual = quality;
  char *backwardQual = NULL;
  char *qual = NULL;
  state_t *state, *currentState;
  int strand, rid;
  bwtint_t nHits = 0;
  unsigned int hitId = 0;
  unsigned int nErrors = states->minErrors[depth];
  size_t nStates = states->nStates[depth][nErrors] = simplifyStates(states->states[depth][nErrors], states->nStates[depth][nErrors]);
  state_t *theseStates = states->states[depth][nErrors];
  writeQname(qname, counts);
  for (size_t i = 0; i < nStates; ++i) {
    nHits += theseStates[i].interval.l - theseStates[i].interval.k + 1;
  }
  if (nHits > parameters->maxNHits) {
    fprintf(outputSamFile, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tNH:i:%lu\tNM:i:%u\n", qname, forwardSeq, forwardQual, nHits, nErrors);
    return;
  }
  if ((nHits == 1) && (nErrors == 0)) {
    pos = bwa_sa2pos(bns, bwt, theseStates[0].interval.k, depth, &strand);
    rid = bns_pos2rid(bns, pos);
    pos = pos - bns->anns[rid].offset + 1;
    printReadUniqueNoError(strand, bns->anns[rid].name, pos, readLength, forwardSeq, quality, counts, outputSamFile);
    return;
  }
  if (nErrors == 0) {
    sprintf(forwardCigar, "%zu=", readLength);
    strcpy(backwardCigar, forwardCigar);
  }
  if ((nHits > 1) || (nErrors >= 40)) {
    mapq = 0;
  }
  else {
    mapq = 40 - nErrors;
  }
  //printf("Read: %s, read length: %zu\n", forwardSeq, readLength);
  for (size_t stateId = 0; stateId < nStates; ++stateId) {
    state = &theseStates[stateId];
    if (nErrors != 0) {
      currentState = state;
      backtraceSize = 0;
      while (currentState->previousState != NULL) {
        //printState(currentState, path->maxDepth);
        char cigar = CIGAR[currentState->trace >> BACKTRACE_OFFSET];
        if ((backtraceSize == 0) || (backtraceCigar[backtraceSize-1] != cigar)) {
          backtraceCigar[backtraceSize] = cigar;
          backtraceLengths[backtraceSize] = 1;
          ++backtraceSize;
        }
        else {
          ++backtraceLengths[backtraceSize-1];
        }
        currentState = currentState->previousState;
      }
      //TODO Optimize this
      forwardCigarSize = backwardCigarSize = 0;
      for (size_t backtraceId = 0; backtraceId < backtraceSize; ++backtraceId) {
        forwardCigarSize  += sprintf(forwardCigar+forwardCigarSize,   "%i%c", backtraceLengths[backtraceId], backtraceCigar[backtraceId]);
        backwardCigarSize += sprintf(backwardCigar+backwardCigarSize, "%i%c", backtraceLengths[backtraceSize-backtraceId-1], backtraceCigar[backtraceSize-backtraceId-1]);
      }
      forwardCigar[forwardCigarSize] = 0;
      backwardCigar[backwardCigarSize] = 0;
    }
    for (bwtint_t j = getStateInterval(state)->k; j <= getStateInterval(state)->l; ++j) {
      flag = (hitId == 0)? 0: CIGAR_SECONDARY_HIT;
      ++hitId;
      pos = bwa_sa2pos(bns, bwt, j, depth, &strand);
      rid = bns_pos2rid(bns, pos);
      pos = pos - bns->anns[rid].offset + 1;
      if (strand == 0) {
        flag |= CIGAR_REVERSE;
        if (backwardSeq  == NULL) backwardSeq  = reverseComplementSequence(forwardSeq, readLength);
        if (backwardQual == NULL) backwardQual = reverseSequence(forwardQual, readLength);
        seq  = backwardSeq;
        qual = backwardQual;
        cigar = backwardCigar;
      }
      else {
        seq  = forwardSeq;
        qual = forwardQual;
        cigar = forwardCigar;
      }
      fprintf(outputSamFile, "%s\t%u\t%s\t%" PRId64 "\t%d\t%s\t*\t0\t0\t%s\t%s\tNH:i:%lu\tHI:i:%u\tIH:i:%lu\tNM:i:%u\n", qname, flag, bns->anns[rid].name, pos, mapq, cigar, seq, qual, nHits, hitId, nHits, nErrors);
    }
  }
  if (backwardSeq)  free(backwardSeq);
  if (backwardQual) free(backwardQual);
}

bool mapWithoutError (states_t *states, size_t depth, unsigned short nt, size_t nErrors) {
  assert(depth > 0);
  state_t *previousState;
  state_t nextState;
  bool mapFound = false;
  //printf("    Mapping %hu without error at depth %zu with %zu errors and %zu states\n", nt, depth, nErrors, states->nStates[depth-1][nErrors]);
  for (size_t stateId = 0; stateId < states->nStates[depth-1][nErrors]; ++stateId) {
    previousState = &states->states[depth-1][nErrors][stateId];
    if (goDownBwt(states, previousState, nt, &nextState)) {
      mapFound = true;
      nextState.trace        |= MATCH;
      nextState.previousState = previousState;
      //printState(&nextState, depth);
      if (! addState(states, depth, nErrors, &nextState)) {
        //printf("      cannot add state\n");
        return false;
      }
    }
  }
  //printf("    found map: %s\n", mapFound ? "true" : "false");
  return mapFound;
}

bool _addError (states_t *states, path_t *path, size_t nErrors, size_t depth, state_t *newState) {
  //printf("    depth %zu, %zu errors, nucleotide %c, %zu states, min errors: %zu\n", depth, nErrors, "ACGT"[path->nucleotides[depth-1]],  states->nStates[depth-1][nErrors-1], states->minErrors[depth-1]);
  //printf("      Path is: ");
  //for (size_t i = 0; i < depth; ++i) putchar("ACGT"[path->nucleotides[i]]);
  //putchar('\n');
  //TODO may be skipped?
  if ((states->maxErrors[depth] != SIZE_MAX) && (states->maxErrors[depth] >= nErrors)) {
    //printf("      first case: %zu/%zu/%zu %zu/%i\n", states->maxErrors[depth], nErrors, SIZE_MAX, states->nStates[depth][nErrors], N_STATES);
    return (states->nStates[depth][nErrors] < N_STATES);
  }
  for (size_t stateId = 0; stateId < states->nStates[depth-1][nErrors-1]; ++stateId) {
    state_t state = states->states[depth-1][nErrors-1][stateId];
    //printState(state, path->maxDepth);
    state.trace         = (state.trace & PREPROCESSED) | INSERTION;
    state.previousState = &states->states[depth-1][nErrors-1][stateId];
    if (! addState(states, depth, nErrors, &state)) {
      //printf("      second case\n");
      return false;
    }
    for (unsigned short nt = 0; nt < N_NUCLEOTIDES; ++nt) {
      if (goDownBwt(states, &state, nt, newState)) {
        //addState(states, depth-1, nErrors, &newState);
        if (nt != path->nucleotides[depth-1]) {
          //printState(newState, path->maxDepth);
          newState->trace        |= MISMATCH | nt;
          newState->previousState = &states->states[depth-1][nErrors-1][stateId];
          if (! addState(states, depth, nErrors, newState)) {
            //printf("      third case\n");
            return false;
          }
        }
      }
    }
  }
  for (size_t stateId = 0; stateId < states->nStates[depth][nErrors-1]; ++stateId) {
    state_t *state = &states->states[depth][nErrors-1][stateId];
    for (unsigned short nt = 0; nt < N_NUCLEOTIDES; ++nt) {
      if (goDownBwt(states, state, nt, newState)) {
        newState->trace        |= DELETION | nt;
        newState->previousState = state;
        if (! addState(states, depth, nErrors, newState)) {
          return false;
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
  //states->nStates[depth][nErrors] = simplifyStates(states->states[depth][nErrors], states->nStates[depth][nErrors]);
  if (states->nStates[depth][nErrors] >= N_STATES) {
    //printf("      fifth case: # states = %zu >= %i after simplification.\n", states->nStates[depth][nErrors], N_STATES);
    return false;
  }
  return true;
}

bool addError (states_t *states, path_t *path) {
  state_t newState;
  size_t nErrors = states->minErrors[path->depth-1] + 1;
  size_t firstDepth;
  bool firstDepthFound = false;
  for (firstDepth = path->depth; (firstDepth > 1) && (! firstDepthFound); --firstDepth) {
    for (size_t nucleotide = 0; (nucleotide < N_NUCLEOTIDES) && (! firstDepthFound); ++nucleotide) {
      firstDepthFound = (states->nStates[firstDepth][nErrors] > 0);
      //if (firstDepthFound) printf("Found first depth at depth = %zu, # errors = %zu, nucleotide = %zu\n", firstDepth, nErrors, nucleotide);
    }
  }
  if (firstDepthFound) firstDepth += 2;
  //printf("  adding error with %zu @ %zu from %zu\n", nErrors, path->depth, firstDepth);
  for (size_t depth = firstDepth; depth <= path->depth; ++depth) {
    if (! _addError(states, path, nErrors, depth, &newState)) {
      //printf("    too many errors\n");
      return false;
    }
  }
  return true;
}

bool shortCutCondition (const states_t *states, const path_t *path) {
  return ((path->depth >= path->shortCut->depth) && (states->minErrors[path->depth] == 0) && (states->nStates[path->depth][0] == 1) && (getStateInterval(&states->states[path->depth][0][0])->k == getStateInterval(&states->states[path->depth][0][0])->l));
}

bool tryShortCut (const tree_t *tree, states_t *states, path_t *path, FILE *outputSamFile) {
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

bool findBestMapping (states_t *states, path_t *path) {
  //printf("  Finding best mapping\n");
  //printf("    Path is: ");
  //for (size_t i = 0; i < path->depth; ++i) putchar("ACGT"[path->nucleotides[i]]);
  //putchar('\n');
  //printStates(states, path->depth);
  //for (size_t i = 0; i < states->nStates[path->depth-1][states->minErrors[path->depth-1]]; ++i) { printState(&states->states[path->depth-1][states->minErrors[path->depth-1]][i], path->maxDepth); }
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

void _map (const tree_t *tree, states_t *states, path_t *path, FILE *outputSamFile) {
  bool mappable = true;
  char *quality;
  while (true) {
    if (! goNextTree(tree, states, path, mappable)) {
      return;
    }
    mappable = findBestMapping(states, path);
    if (mappable) {
      if ((quality = findQuality(&tree->qualities, path->cellIds[path->depth])) != NULL) {
        printRead(states, path, quality, tree->cells[path->cellIds[path->depth]].counts, outputSamFile);
      }
      if (shortCutCondition(states, path)) {
        if (! tryShortCut(tree, states, path, outputSamFile)) {
          //printf("Short cut with negative exit\n");
          return;
        }
      }
    }
  }
}

FILE *printSamHeader() {
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

void map (const tree_t *tree, FILE *outputSamFile) {
  states_t *states = initializeStates(tree->depth);
  path_t   *path   = initializePath(tree->depth);
  //addState(&states, 0, 0, firstState);
  //goDownTree(tree, &path);
  _map(tree, states, path, outputSamFile);
  freeStates(states);
  freePath(path);
}

int main(int argc, char const ** argv) {
  int returnCode = 0;
  parameters_t param;
  stats_t stat;
  tree_t tree;
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
  FILE *outputSamFile = printSamHeader();
  if (outputSamFile == NULL) {
    printf("Error!  Cannot write to output SAM file '%s'.\nExiting.\n", param.outputSamFileName);
  }
  map(&tree, outputSamFile);
  freeTree(&tree);
  bwa_idx_destroy(idx);
  puts("... done.");
  printStats();
  return 0;
}
