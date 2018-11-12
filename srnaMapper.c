// ==========================================================================
//                              srnaMapper
// ==========================================================================
// Copyright (C) 2018 Matthias Zytnicki, INRA
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
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
// Parse parameters
// ==========================================================================

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>

#include "Libs/bwa/bwt.h"
#include "Libs/bwa/bwa.h"
#include "Libs/bwa/bwase.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

#define N_NUCLEOTIDES    4
#define NUCLEOTIDES_BITS 2
#define TRIPLET          3
#define N_TRIPLETS      64
#define TRIPLET_MASK    15

#define MAX_HITS        20

#define INIT_N_CELLS 0x1000000

#define NO_DATA        0

typedef unsigned long count_t;

const char DNA5_TO_CHAR [5] = { 'A', 'C', 'G', 'T', 'N' };
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

void trimSequence (size_t l, char *s) {
  s[l-1] = 0;
}


typedef struct {
    char *readsFileNames[255];
    char *genomeFileName;
    char *outputReadsFileName;
    unsigned int nReadsFiles;
    unsigned int lowComplexityThreshold;
} parameters_t;

parameters_t *parameters;

void printUsage () {
  puts("srnaCollapser [-h] -r reads -g genome [-o filename] [-f filter]");
}

int parseCommandLine (int argc, char const **argv) {
  char *endptr;
  parameters->outputReadsFileName    = NULL;
  parameters->nReadsFiles            = 0;
  parameters->lowComplexityThreshold = 6;
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
      parameters->outputReadsFileName = strdup(argv[i]);
    }
    else if (strcmp(argv[i], "-f") == 0) {
      ++i;
      parameters->lowComplexityThreshold = strtol(argv[i], &endptr, 10);
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
  return EXIT_SUCCESS;
}

typedef struct {
  uint64_t children [N_NUCLEOTIDES];
  char    *quality;
  count_t *counts;
} cell_t;

void createCell (cell_t *cell) {
  for (unsigned short i = 0; i < N_NUCLEOTIDES; ++i) {
    cell->children[i] = NO_DATA;
  }
  cell->counts = (count_t *) calloc(parameters->nReadsFiles, sizeof(count_t));
  cell->quality = NULL;
}

void freeCell (cell_t *cell) {
  if (cell->quality) free(cell->quality);
  free(cell->counts);
}

void setQuality (cell_t *cell, size_t l, char *quality, unsigned int fileId) {
  ++cell->counts[fileId];
  if (cell->quality == NULL) {
    cell->quality = strdup(quality);
  }
  else {
    for (size_t i = 0; i < l; ++i) {
      cell->quality[i] = MAX(cell->quality[i], quality[i]);
    }
  }
}

typedef struct {
  size_t   depth;
  cell_t  *cells;
  uint64_t nCells;
  uint64_t nAllocated;
} tree_t;

void createTree (tree_t *tree) {
  tree->nAllocated = INIT_N_CELLS;
  tree->cells = (cell_t *) malloc(tree->nAllocated * sizeof(cell_t));
  tree->depth = 0;
  tree->nCells = 0;
  createCell(&tree->cells[0]);
}

void freeTree (tree_t *tree) {
  for (uint64_t i = 0; i < tree->nCells; ++i) {
    freeCell(&tree->cells[i]);
  }
  free(tree->cells);
}

uint64_t addCell (tree_t *tree) {
  if (tree->nCells == tree->nAllocated) {
    tree->nAllocated *= 2;
    if ((tree->cells = (cell_t *) realloc(tree->cells, tree->nAllocated * sizeof(cell_t))) == NULL) {
      printf("Cannot allocate memory for tree of size %lu.\nExiting.\n", tree->nAllocated);
      exit(EXIT_FAILURE);
    }
  }
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

void addSequence (tree_t *tree, size_t l, char *sequence, char *quality, unsigned int fileId) {
  uint64_t cellId = 0;
  for (int sequenceId = l-1; sequenceId >= 0; --sequenceId) {
    cellId = goDown(tree, cellId, CHAR_TO_DNA5[(int) sequence[sequenceId]]);
  }
  setQuality(&tree->cells[cellId], l, quality, fileId);
  tree->depth = MAX(tree->depth, l);
}

void _printTree (const tree_t *tree, FILE *outFile, uint64_t *readId, char *read, size_t readPos, uint64_t cellId, unsigned short triplet, count_t *tripletCount) {
  uint64_t nextCellId;
  unsigned short nextTriplet;
  cell_t *cell = &tree->cells[cellId];
  //read[readPos] = 0; printf("Base triplet: %d, size: %zu, read: %s\n", triplet, readPos, read);
  //read[readPos] = 0; printf("Current state: %zu, %s, %zu, %lu\n", *readId, read, readPos, cellId); fflush(stdout);
  if (cell->quality) {
    //printf("\tGot quality\n"); fflush(stdout);
    ++(*readId);
    read[readPos] = 0;
    fprintf(outFile, "@read%"PRIu64"_x", *readId);
    for (unsigned int fileId = 0; fileId < parameters->nReadsFiles; ++fileId) {
      fprintf(outFile, "_%lu", cell->counts[fileId]);
    }
    fprintf(outFile, "\n%s\n+\n%s\n", read, cell->quality);
  }
  triplet &= TRIPLET_MASK;
  triplet <<= NUCLEOTIDES_BITS;
  //printf("\tNew triplet: %d\n", triplet);
  for (unsigned short i = 0; i < N_NUCLEOTIDES; ++i) {
    //printf("\tTrying child %d\n", i); fflush(stdout);
    nextCellId = cell->children[i];
    if (nextCellId != NO_DATA) {
      nextTriplet = triplet | i;
      //printf("\tNext triplet (%d): %d\n", i, nextTriplet);
      if (readPos >= TRIPLET-1) tripletCount[nextTriplet] += 1;
      if (tripletCount[nextTriplet] <= parameters->lowComplexityThreshold) {
        read[readPos] = DNA5_TO_CHAR[i];
        _printTree(tree, outFile, readId, read, readPos+1, nextCellId, nextTriplet, tripletCount);
      }
      if (readPos >= TRIPLET-1) tripletCount[nextTriplet] -= 1;
    }
  }
}

int printTree (char *fileName, const tree_t *tree) {
  FILE *outFile = fopen(fileName, "w");
  if (outFile == NULL) return EXIT_FAILURE;
  uint64_t readId = 0;
  char *read = (char *) malloc((tree->depth+1) * sizeof(char));
  count_t tripletCount [N_TRIPLETS];
  for (unsigned int i = 0; i < N_TRIPLETS; ++i) tripletCount[i] = 0;
  _printTree(tree, outFile, &readId, read, 0, 0, 0, tripletCount);
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
    getline(&sequence, &len, inFile);
    getline(&line, &len, inFile);
    nRead = getline(&quality, &len, inFile);
    trimSequence(nRead, sequence);
    trimSequence(nRead, quality);
    addSequence(tree, nRead-1, sequence, quality, fileId);
  }
  free(line);
  free(sequence);
  free(quality);
  fclose(inFile);
  return EXIT_SUCCESS;
}

bwaidx_t *loadGenomeFile (char *indexName) {
  bwaidx_t *idx = bwa_idx_load(indexName, BWA_IDX_ALL);
  if (idx == NULL) {
    fprintf(stderr, "Index load failed.\n");
  }
  return idx;
}

void _map (const bwt_t *bwt, const bntseq_t *bns, bwtint_t k, bwtint_t l, const tree_t *tree, uint64_t cellId, size_t readPos, char *read, unsigned short triplet, count_t *tripletCount) {
  bwtint_t nk, nl, ok, ol;
  cell_t *cell = &tree->cells[cellId];
  int64_t pos;
  uint64_t nextCellId;
  int strand, rid;
  unsigned short nextTriplet;
  if (cell->quality) {
    printf("Mapped '%s'\n", read+tree->depth-readPos);
    if (l - k + 1 > MAX_HITS) {
      printf("\t%lu hits\n", l - k + 1);
    }
    else {
      for (bwtint_t i = k; i <= l; ++i) {
        pos = bwa_sa2pos(bns, bwt, i, readPos, &strand);
        rid = bns_pos2rid(bns, pos);
        pos = pos - bns->anns[rid].offset;
        printf("\t%d %" PRId64 " (%c)\n", rid, pos, "-+"[strand]);
      }
    }
  }
  triplet &= TRIPLET_MASK;
  triplet <<= NUCLEOTIDES_BITS;
  for (unsigned short i = 0; i < N_NUCLEOTIDES; ++i) {
    nextCellId = cell->children[i];
    if (nextCellId != NO_DATA) {
      nextTriplet = triplet | i;
      if (readPos >= TRIPLET-1) tripletCount[nextTriplet] += 1;
      if (tripletCount[nextTriplet] <= parameters->lowComplexityThreshold) {
        bwt_2occ(bwt, k-1, l, i, &ok, &ol);
        nk = bwt->L2[i] + ok + 1;
        nl = bwt->L2[i] + ol;
        if (nk <= nl) {
          read[tree->depth-readPos-1] = DNA5_TO_CHAR[i];
          _map(bwt, bns, nk, nl, tree, nextCellId, readPos+1, read, nextTriplet, tripletCount);
        }
      }
      if (readPos >= TRIPLET-1) tripletCount[nextTriplet] -= 1;
    }
  }
}

void map (const bwt_t *bwt, const bntseq_t *bns, const tree_t *tree) {
  char *read = (char *) malloc((tree->depth+1) * sizeof(char));
  read[tree->depth] = 0;
  count_t tripletCount [N_TRIPLETS];
  for (unsigned int i = 0; i < N_TRIPLETS; ++i) tripletCount[i] = 0;
  _map(bwt, bns, 0, bwt->seq_len, tree, 0, 0, read, 0, tripletCount);
  free(read);
}

int main(int argc, char const ** argv) {
  int returnCode = 0;
  parameters_t param;
  tree_t tree;
  bwaidx_t *idx = NULL;
  parameters = &param;
  returnCode = parseCommandLine(argc, argv);
  if (returnCode != EXIT_SUCCESS) return returnCode;
  createTree(&tree);
  for (unsigned int fileId = 0; fileId < parameters->nReadsFiles; ++fileId) {
    printf("Reading reads file %s...\n", parameters->readsFileNames[fileId]);
    returnCode = readReadsFile(parameters->readsFileNames[fileId], &tree, fileId);
    if (returnCode != EXIT_SUCCESS) return returnCode;
    puts("... done.");
  }
  printf("Depth of tree: %zu\n", tree.depth);
  if (parameters->outputReadsFileName != NULL) {
    puts("Printing tree...");
    returnCode = printTree(parameters->outputReadsFileName, &tree);
    if (returnCode != EXIT_SUCCESS) return returnCode;
    puts("... done.");
  }
  idx = loadGenomeFile(parameters->genomeFileName);
  if (idx == NULL) return EXIT_FAILURE;
  puts("... done.");
  map(idx->bwt, idx->bns, &tree);
  freeTree(&tree);
  bwa_idx_destroy(idx);
  return 0;
}
