#ifndef BWT_H
#define BWT_H

#include "constants.h"
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
  // TODO: Does that really work?  Maybe a padding problem?
  // return (memcmp(i1, i2, sizeof(bwtinterval_t)) == 0);
  return ((i1->k == i2->k) && (i1->l == i2->l));
}


size_t getBwtHash (bwtinterval_t interval) {
  //TODO: Compute a "good" hash function here
  // return ((interval.k * 98317 ^ interval.l * 49157) % BWT_BUFFER_SIZE);
  return (interval.k + interval.l);
}

bool goDownBwt (state_t *previousState, unsigned short nucleotide, bwtinterval_t *newInterval) {
  bwt_2occ(bwt, previousState->interval.k-1, previousState->interval.l, nucleotide, &newInterval->k, &newInterval->l);
  newInterval->k = bwt->L2[nucleotide] + newInterval->k + 1;
  newInterval->l = bwt->L2[nucleotide] + newInterval->l;
  //if (newState->interval.k <= newState->interval.l) printf("      ok!\n");
  return (newInterval->k <= newInterval->l);
}

#endif
