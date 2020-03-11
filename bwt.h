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
  assert(i1->l != 0);
  //if (i2->l == 0) return false;
  //if ((i1->l == 0) || (i2->l == 0)) return false;
  return ((i1->k == i2->k) && (i1->l == i2->l));
}


typedef struct {
  bwtinterval_t interval;
  bwtinterval_t nextIntervals[N_NUCLEOTIDES];
} bwt_buffer_data_t;

typedef struct {
  bwt_buffer_data_t intervals[BWT_BUFFER_SIZE];
} bwt_buffer_t;

void createBwtBuffer (bwt_buffer_t *bwtBuffer) {
  for (size_t bufferId = 0; bufferId < BWT_BUFFER_SIZE; ++bufferId) {
    bwtBuffer->intervals[bufferId].interval.l = 0;
    for (unsigned short nt = 0; nt < N_NUCLEOTIDES; ++nt) {
      bwtBuffer->intervals[bufferId].nextIntervals[nt].l = 0;
    }
  }
}

size_t getBwtBufferIndex (bwtinterval_t interval) {
  //TODO: Compute a "good" hash function here
  //return (interval.k * 106039 + interval.l) & BWT_BUFFER_MASK;
  return ((interval.k + interval.l) % BWT_BUFFER_SIZE);
  //return ((interval.k + interval.l) % BWT_BUFFER_SIZE);
}

void addToBwtBuffer (bwt_buffer_t *bwtBuffer, bwtinterval_t previousInterval, unsigned short nucleotide, bwtinterval_t nextInterval) {
  size_t bufferId = getBwtBufferIndex (previousInterval);
  bwt_buffer_data_t *bufferData = &bwtBuffer->intervals[bufferId];
  if (! compareBwtIntervals(&previousInterval, &bufferData->interval)) {
    bufferData->interval = previousInterval;    
    for (unsigned short nt = 0; nt < N_NUCLEOTIDES; ++nt) {
      bufferData->nextIntervals[nt].l = 0;
    }
  }
  bufferData->nextIntervals[nucleotide] = nextInterval;
}

bwtinterval_t *findInBwtBuffer (bwt_buffer_t *bwtBuffer, bwtinterval_t previousInterval, unsigned short nucleotide) {
  size_t bufferId = getBwtBufferIndex (previousInterval);
  bwt_buffer_data_t *bufferData = &bwtBuffer->intervals[bufferId];
  ++stats->nBufferCalls;
  if (compareBwtIntervals(&previousInterval, &bufferData->interval)) {
    if (bufferData->nextIntervals[nucleotide].l != 0) {
      ++stats->nBufferCallSucesses;
      //printf("Buffered!\n");
      return &bufferData->nextIntervals[nucleotide];
    }
  }
  return NULL;
}

/*
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
*/

bool goDownBwt (bwt_buffer_t *bwtBuffer, state_t *previousState, unsigned short nucleotide, bwtinterval_t *newInterval) {
  bwtinterval_t *bufferedInterval = findInBwtBuffer(bwtBuffer, previousState->interval, nucleotide);
  if (bufferedInterval != NULL) {
    *newInterval = *bufferedInterval;
    return (newInterval->k <= newInterval->l);
  }
  ++stats->nDown;
  bwt_2occ(bwt, previousState->interval.k-1, previousState->interval.l, nucleotide, &newInterval->k, &newInterval->l);
  newInterval->k = bwt->L2[nucleotide] + newInterval->k + 1;
  newInterval->l = bwt->L2[nucleotide] + newInterval->l;
  //if (newState->interval.k <= newState->interval.l) printf("      ok!\n");
  addToBwtBuffer(bwtBuffer, previousState->interval, nucleotide, *newInterval);
  return (newInterval->k <= newInterval->l);
}

/*
void goDownBwt3Nt (bwt_buffer_t *bwtBuffer, state_t *previousState, unsigned short nucleotide, bwtinterval_t *newIntervals) {
  for (unsigned short nt = 0; nt < N_NUCLEOTIDES; ++nt) {
    ++stats->nDown;
    if (nt != nucleotide) {
      bwt_2occ(bwt, previousState->interval.k-1, previousState->interval.l, nucleotide, &newIntervals[nt].k, &newIntervals[nt].l);
      newIntervals[nt].k = bwt->L2[nt] + newIntervals[nt].k + 1;
      newIntervals[nt].l = bwt->L2[nt] + newIntervals[nt].l;
    }
  }
}

void goDownBwt4Nt (bwt_buffer_t *bwtBuffer, state_t *previousState, bwtinterval_t *newIntervals) {
  for (unsigned short nt = 0; nt < N_NUCLEOTIDES; ++nt) {
    ++stats->nDown;
    bwt_2occ(bwt, previousState->interval.k-1, previousState->interval.l, nt, &newIntervals[nt].k, &newIntervals[nt].l);
    newIntervals[nt].k = bwt->L2[nt] + newIntervals[nt].k + 1;
    newIntervals[nt].l = bwt->L2[nt] + newIntervals[nt].l;
  }
}
*/


#endif
