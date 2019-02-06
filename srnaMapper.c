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

#include <stdbool.h>
#include <assert.h>
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

#define MAX_PATHS       0x10
#define MAX_HITS        20
#define N_STATES        0x10000

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
    size_t       maxNErrors;
    unsigned int lowComplexityThreshold;
} parameters_t;

parameters_t *parameters;

void printUsage () {
  puts("srnaCollapser [-h] -r reads -g genome [-o filename] [-f filter] [-e #errors]");
}

int parseCommandLine (int argc, char const **argv) {
  char *endptr;
  parameters->outputReadsFileName    = NULL;
  parameters->nReadsFiles            = 0;
  parameters->maxNErrors             = 2;
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
    else if (strcmp(argv[i], "-e") == 0) {
      ++i;
      parameters->maxNErrors = strtol(argv[i], &endptr, 10);
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
  //printf("%u %zu %p\n", parameters->nReadsFiles, sizeof(count_t), cell);
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
    assert(strlen(cell->quality) == strlen(quality));
    assert(strlen(cell->quality) == l);
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
  tree->nCells = 1;
  createCell(&tree->cells[0]);
}

void freeTree (tree_t *tree) {
  tree->nAllocated = 0;
  tree->depth = 0;
  tree->nCells = 0;
  for (uint64_t i = 0; i < tree->nCells; ++i) {
    freeCell(&tree->cells[i]);
  }
  free(tree->cells);
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

void addSequence (tree_t *tree, size_t l, char *sequence, char *quality, unsigned int fileId) {
  uint64_t cellId = 0;
  assert(strlen(sequence) == strlen(quality));
  assert(strlen(quality) == l);
  for (int sequenceId = l-1; sequenceId >= 0; --sequenceId) {
    cellId = goDown(tree, cellId, CHAR_TO_DNA5[(int) sequence[sequenceId]]);
  }
  setQuality(&tree->cells[cellId], l, quality, fileId);
  tree->depth = MAX(tree->depth, l);
}

void _printTree (const tree_t *tree, FILE *outFile, uint64_t *readId, char *read, size_t readPos, uint64_t cellId) {
  uint64_t nextCellId;
  cell_t *cell = &tree->cells[cellId];
  //read[readPos] = 0; printf("Base triplet: %d, size: %zu, read: %s\n", triplet, readPos, read);
  //read[readPos] = 0; printf("Current state: %zu, %s, %zu, %lu\n", *readId, read, readPos, cellId); fflush(stdout);
  if (cell->quality) {
    //printf("\tGot quality\n"); fflush(stdout);
    ++(*readId);
    fprintf(outFile, "@read%"PRIu64"_x", *readId);
    for (unsigned int fileId = 0; fileId < parameters->nReadsFiles; ++fileId) {
      fprintf(outFile, "_%lu", cell->counts[fileId]);
    }
    fprintf(outFile, "\n%s\n+\n%s\n", read+tree->depth-readPos, cell->quality);
    assert(strlen(read+tree->depth-readPos) == strlen(cell->quality));
  }
  //printf("\tNew triplet: %d\n", triplet);
  for (unsigned short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    //printf("\tTrying child %d\n", nucleotide); fflush(stdout);
    nextCellId = cell->children[nucleotide];
    if (nextCellId != NO_DATA) {
      //printf("\tNext triplet (%d): %d\n", i, nextTriplet);
      read[tree->depth-readPos-1] = DNA5_TO_CHAR[nucleotide];
      _printTree(tree, outFile, readId, read, readPos+1, nextCellId);
    }
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
    getline(&sequence, &len, inFile);
    getline(&line, &len, inFile);
    nRead = getline(&quality, &len, inFile);
    assert(strlen(sequence) == strlen(quality));
    assert(strlen(sequence) == (unsigned long) nRead);
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

bool _filterTree (const tree_t *tree, size_t readPos, uint64_t cellId, unsigned short triplet, count_t *tripletCount) {
  uint64_t nextCellId;
  unsigned short nextTriplet;
  bool foundRead = false;
  cell_t *cell = &tree->cells[cellId];
  triplet &= TRIPLET_MASK;
  triplet <<= NUCLEOTIDES_BITS;
  if (cell->quality) {
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
        if (! _filterTree(tree, readPos+1, nextCellId, nextTriplet, tripletCount)) {
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
} state_t;

typedef struct {
  state_t ***states;
  size_t **nStates;
  size_t *nStatesPerPosition;
  size_t *minErrors;
  size_t treeSize;
} states_t;

bool compareStates (state_t *state1, state_t *state2) {
  return ((state1->k == state2->k) && (state1->l == state2->l));
}

void addState(states_t *states, size_t depth, size_t nErrors, state_t *state) {
  //printf("\t\t\tAdding one state (%" PRIu64 ", %" PRIu64 ") at (%zu, %zu), %zu/%d occupied\n", state->k, state->l, depth, nErrors, states->nStates[depth][nErrors], N_STATES);
  for (size_t stateId = 0; stateId < states->nStates[depth][nErrors]; ++stateId) {
    if (compareStates(state, &states->states[depth][nErrors][stateId])) {
      return;
    }
  }
  //printf("\t\t\t\tadding it\n");
  assert(states->nStates[depth][nErrors] < N_STATES-1);
  states->states[depth][nErrors][states->nStates[depth][nErrors]] = *state;
  ++states->nStates[depth][nErrors];
  ++states->nStatesPerPosition[nErrors];
  if (nErrors < states->minErrors[depth]) states->minErrors[depth] = nErrors;
}

state_t *initializeStates(states_t *states, size_t treeSize, bwtint_t genomeSize) {
  states->treeSize           = treeSize;
  states->states             = (state_t ***) malloc(treeSize * sizeof(states_t **));
  states->nStates            = (size_t **)   malloc(treeSize * sizeof(size_t *));
  states->nStatesPerPosition = (size_t *)    calloc(treeSize,  sizeof(size_t));
  states->minErrors          = (size_t *)    malloc(treeSize * sizeof(size_t));
  for (size_t i = 0; i < treeSize; ++i) {
    states->states[i]    = (state_t **) malloc((parameters->maxNErrors+1) * sizeof(states_t *));
    states->nStates[i]   = (size_t *)   calloc((parameters->maxNErrors+1),  sizeof(size_t));
    states->minErrors[i] = SIZE_MAX;
    for (size_t j = 0; j <= parameters->maxNErrors; ++j) {
      states->states[i][j] = (state_t *) malloc(N_STATES * sizeof(states_t));
    }
  }
  state_t baseState = { 0, genomeSize };
  addState(states, 0, 0, &baseState);
  return &(states->states[0][0][0]);
}

/*
bool getStates(states_t *states, size_t readPos, size_t *nCurrentStates, state_t **currentStates, short int *nErrors) {
  if (states->nStatesPerPosition[readPos] == 0) {
    return false;
  }
  for (size_t j = 0; j <= maxNErrors; ++j) {
    if (states->nStates[readPos][j] > 0) {
      *nErrors        = j;
      *nCurrentStates = nStates[readPos][j];
      *currentStates  = states[readPos][j];
      return true;
    }
  }
  assert(false);
  return false;
}
*/

void backtrackStates(states_t *states, size_t level) {
  //printf("\t\t\tBacktracking to level %zu\n", level);
  for (size_t i = level+1; i < states->treeSize; ++i) {
    if (states->nStatesPerPosition[i] == 0) {
      return;
    }
    states->nStatesPerPosition[i] = 0;
    states->minErrors[i]          = SIZE_MAX;
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
}

typedef struct {
  unsigned short *nucleotides;
  uint64_t *cellIds;
  size_t maxDepth;
  size_t depth;
  char *read;
  size_t readPos;
} path_t;

void initializePath (path_t *path, size_t maxDepth) {
  path->nucleotides    = malloc(maxDepth * sizeof(unsigned short));
  path->cellIds        = malloc(maxDepth * sizeof(uint64_t));
  path->read           = malloc((maxDepth+1) * sizeof(char));
  path->maxDepth       = maxDepth;
  path->depth          = 0;
  path->cellIds[0]     = 0;
  path->read[maxDepth] = 0;
  path->readPos        = maxDepth;
}

void freePath (path_t *path) {
  free(path->nucleotides);
  free(path->cellIds);
}

bool gotoNextNucleotide (const tree_t *tree, path_t *path) {
  uint64_t cellId;
  //printf("  Goto next nt from height %zu with read '%s' ", path->depth, path->read+path->readPos);
  cell_t *cell = &tree->cells[path->cellIds[path->depth]];
  for (unsigned short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    if ((cellId = cell->children[nucleotide]) != NO_DATA) {
      path->nucleotides[path->depth] = nucleotide;
      ++path->depth;
      --path->readPos;
      path->cellIds[path->depth] = cellId;
      path->read[path->readPos] = DNA5_TO_CHAR[nucleotide];
      //printf("to %zu with read '%s'\n", path->depth, path->read+path->readPos);
      return true;
    }
  }
  //printf("to nothing\n");
  return false;
}

bool gotoNextRead (const tree_t *tree, path_t *path) {
  uint64_t cellId;
  cell_t *cell;
  //printf("  Goto next read from depth %zu with read '%s'\n", path->depth, path->read+path->readPos);
  while (path->depth > 0) {
    --path->depth;
    ++path->readPos;
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

bool mapOneNucleotide (const bwt_t *bwt, states_t *states, unsigned short nucleotide, size_t depth) {
  //TODO probably some paths are missing here
  bwtint_t nk, nl;
  state_t *previousState;
  state_t nextState;
  bool newState = false;
  //printf("\tMapping 1 nt at depth %zu and nt %hu, with at most %zu errors\n", depth, nucleotide, parameters->maxNErrors);
  for (size_t nErrors = 0; nErrors <= parameters->maxNErrors; ++nErrors) {
    //printf("\t\t# errors: %zu, got %zu states\n", nErrors, states->nStates[depth][nErrors]);
    for (size_t stateId = 0; stateId < states->nStates[depth][nErrors]; ++stateId) {
      previousState = &states->states[depth][nErrors][stateId];
      if (nErrors == parameters->maxNErrors) {
        bwt_2occ(bwt, previousState->k-1, previousState->l, nucleotide, &nk, &nl);
        nk = bwt->L2[nucleotide] + nk + 1;
        nl = bwt->L2[nucleotide] + nl;
        if (nk <= nl) {
          nextState = (state_t) { nk, nl };
          addState(states, depth+1, nErrors, &nextState);
          newState = true;
        }
      }
      else {
        addState(states, depth+1, nErrors+1, previousState);
        newState = true;
        for (unsigned short otherNucleotide = 0; otherNucleotide < N_NUCLEOTIDES; ++otherNucleotide) {
          bwt_2occ(bwt, previousState->k-1, previousState->l, otherNucleotide, &nk, &nl);
          nk = bwt->L2[otherNucleotide] + nk + 1;
          nl = bwt->L2[otherNucleotide] + nl;
          //printf("\t\t\tChecking nt %d between range %lu-%lu, got %lu-%lu\n", otherNucleotide, previousState->k-1, previousState->l, nk, nl);
          if (nk <= nl) {
            nextState = (state_t) { nk, nl };
            addState(states, depth+1, nErrors + ((nucleotide==otherNucleotide)? 0: 1), &nextState);
          }
        }
      }
    }
  }
  //printf("\tOver with mapping 1 and returns %s\n", newState ? "true" : "false");
  return newState;
}

void printRead (const bwt_t *bwt, const bntseq_t *bns, size_t nStates, state_t *states, path_t *path) {
  int64_t pos;
  int strand, rid;
  bwtint_t nHits = 0;
  printf("Mapped '%s'\n", path->read + path->readPos);
  for (size_t i = 0; i < nStates; ++i) {
    nHits += states[i].l - states[i].k + 1;
  }
  if (nHits > MAX_HITS) {
    printf("\t%lu hits\n", nHits);
    return;
  }
  for (size_t i = 0; i < nStates; ++i) {
    for (bwtint_t j = states[i].k; j <= states[i].l; ++j) {
      pos = bwa_sa2pos(bns, bwt, j, path->depth, &strand);
      rid = bns_pos2rid(bns, pos);
      pos = pos - bns->anns[rid].offset;
      printf("\t%d %" PRId64 " (%c)\n", rid, pos, "-+"[strand]);
    }
  }
}

void _map (const bwt_t *bwt, const bntseq_t *bns, const tree_t *tree, states_t *states, path_t *path) {
  assert(path->depth > 0);
  while (true) {
    //printf("Starting 'map' at height %zu with read '%s' (min errors: %zu)\n", path->depth, path->read+path->readPos, states->minErrors[path->depth-1]);
    if (mapOneNucleotide(bwt, states, path->nucleotides[path->depth-1], path->depth-1)) {
      //printf("\tquality is %p\n", tree->cells[path->cellIds[path->depth]].quality);
      if (tree->cells[path->cellIds[path->depth]].quality) {
        printRead(bwt, bns, states->nStates[path->depth][states->minErrors[path->depth]], states->states[path->depth][states->minErrors[path->depth]], path);
      }
      //printf("\tGoing to next nt (depth is %zu)\n", path->depth);
      if (! gotoNextNucleotide(tree, path)) {
        //printf("\tThen to next read (depth is %zu)\n", path->depth);
        if (! gotoNextRead(tree, path)) {
          return;
        }
      }
    }
    else {
      //printf("\tGoing to next read\n");
      if (! gotoNextRead(tree, path)) {
        return;
      }
    }
    backtrackStates(states, path->depth);
  }
}

void map (const bwt_t *bwt, const bntseq_t *bns, const tree_t *tree) {
  states_t states;
  path_t   path;
  state_t *firstState = initializeStates(&states, tree->depth, bwt->seq_len);
  initializePath(&path, tree->depth);
  addState(&states, 0, 0, firstState);
  gotoNextNucleotide(tree, &path);
  _map(bwt, bns, tree, &states, &path);
  freeStates(&states);
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
  puts("Filtering tree...");
  filterTree(&tree);
  printf("... done.\n");
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
  puts("Mapping...");
  map(idx->bwt, idx->bns, &tree);
  freeTree(&tree);
  bwa_idx_destroy(idx);
  puts("... done.");
  return 0;
}
