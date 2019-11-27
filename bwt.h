#ifndef BWT_H
#define BWT_H

#include "constants.h"
#include "stats.h"
#include "state.h"

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

bool compareBwtIntervals (bwtinterval_t *i1, bwtinterval_t *i2) {
  return ((i1->k == i2->k) && (i1->l == i2->l));
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

#endif
