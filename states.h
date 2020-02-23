#ifndef STATES_H
#define STATES_H

#include "constants.h"
#include "parameters.h"
#include "state.h"
#include "bwt.h"
#include "sw.h"

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
  state_t     **states;
  size_t      **firstState;
  size_t      **nStates;
  size_t       *allocatedNStates;
  size_t       *nStatesPerPosition;
  size_t       *minErrors;
  size_t       *maxErrors;
  size_t       depth;
  sw_t         *sw;
  bwt_buffer_t *bwtBuffer;
} states_t;

void printStates (const states_t *states, size_t depth) {
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
  for (size_t i = 0; i <= depth; ++i) {
    printf("\t\t\t\tdepth %zu:", i);
    for (size_t nErrors = 0; nErrors <= parameters->maxNErrors; ++nErrors) {
      printf("\t%zu errors: %zu (%zu) cells, starting @ %zu", nErrors, states->nStates[i][nErrors], states->firstState[i+1][nErrors] - states->firstState[i][nErrors], states->firstState[i][nErrors]);
      //assert(states->firstState[i+1][nErrors] - states->firstState[i][nErrors] == states->nStates[i][nErrors]);
    }
    printf("\n");
  }
}

size_t getNStates(const states_t *states, size_t depth, size_t nErrors) {
  return states->nStates[depth][nErrors];
}

state_t *getState(const states_t *states, size_t depth, size_t nErrors, size_t i) {
  //printf("Getting state #%zu at depth %zu with %zu errors: %zu\n", i, depth, nErrors, states->firstState[depth][nErrors]+i);
  //printf("\tfirst state: %zu, next first state: %zu\n", states->firstState[depth][nErrors], states->firstState[depth+1][nErrors]);
  //printStates(states, depth);
  assert(i < getNStates(states, depth, nErrors));
  //assert((depth == states->depth) || (states->firstState[depth+1][nErrors] == 0) || (states->firstState[depth+1][nErrors] >= states->firstState[depth][nErrors]));
  return &states->states[nErrors][states->firstState[depth][nErrors]+i];
}

void simplifyStates (states_t *states, size_t depth, size_t nErrors) {
  return;
  assert(depth < states->depth);
  assert(nErrors <= states->maxErrors[depth]);
  size_t previousNStates = states->nStates[depth][nErrors];
  size_t nextNStates = 0;
  state_t *theseStates = states->states[nErrors] + states->firstState[depth][nErrors];
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

/*
void heapifyStates (state_t *states, size_t nStates) {
  size_t currentId = nStates-1;
  size_t parentId  = getParentStateId(currentId);
  while ((currentId != 0) && (stateComp(&states[parentId], &states[currentId]) > 0)) { 
    swapStates(&states[currentId], &states[parentId]); 
    currentId = parentId;
  }
}
*/

/*
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
*/

state_t *allocateNewState (states_t *states, size_t depth, size_t nErrors) {
  assert(depth <= states->depth);
  //printf("Allocating state @ depth %zu, %zu errors\n", depth, nErrors);
  //printStates(states, depth+1);
  size_t nStatesPerError;
  /*
  if (isStateInserted(states->states[depth][nErrors], states->nStates[depth][nErrors], state)) {
    return NULL;
  }
  */
  //states->states[depth][nErrors][states->nStates[depth][nErrors]] = *state;
  //printf("# states so far: %zu\n", states->nStates[depth][nErrors]);
  if (states->nStates[depth][nErrors] == 0) {
    // TODO: BAD BUG FIX!  SW should be complete!
    if (depth == 0) {
      states->firstState[depth][nErrors] = 0;
    }
    else {
      size_t previousDepth;
      for (previousDepth = depth-1; (previousDepth > 0) && (states->nStates[previousDepth][nErrors] == 0); --previousDepth) ;
      states->firstState[depth][nErrors] = states->firstState[previousDepth][nErrors] + states->nStates[previousDepth][nErrors];
    }
    //printf("First state of %zu/%zu is now %zu\n", depth, nErrors, states->firstState[depth][nErrors]);
    //states->firstState[depth][nErrors] = states->nStatesPerError[nErrors];
  }
  //printf("first state: %zu\n", states->firstState[depth][nErrors]);
  ++states->nStates[depth][nErrors];
  ++states->nStatesPerPosition[depth];
  nStatesPerError = states->firstState[depth][nErrors] + states->nStates[depth][nErrors];
  stats->maxNStates = MAX(stats->maxNStates, nStatesPerError);
  if (nStatesPerError == states->allocatedNStates[nErrors]-1) {
    states->allocatedNStates[nErrors] *= 2;
    states->states[nErrors] = (state_t *) realloc(states->states[nErrors], states->allocatedNStates[nErrors] * sizeof(state_t));
    if (states->states[nErrors] == NULL) {
      fprintf(stderr, "Memory error: Cannot allocate an array of States for errors %zu of size %zu * %zu\n", nErrors, states->allocatedNStates[nErrors], sizeof(state_t));
      exit(EXIT_FAILURE);
    }
    //printf("Reallocation #states with %zu errors: now %zu\n", nErrors, states->allocatedNStates[nErrors]);
    //printStates(states, depth);
    //return NULL;
  }
  //heapifyStates(states->states[depth][nErrors], states->nStates[depth][nErrors]);
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
  //printf("Add states finishes with index (count 1) %zu\n", nStatesPerError - 1);
  //printStates(states, depth+1);
  //printf("\tSetting at position %zu %p (%zu/%zu)\n", nStatesPerError - 1, &states->states[nErrors][nStatesPerError - 1], states->firstState[depth][nErrors], states->nStates[depth][nErrors]);
  return &states->states[nErrors][nStatesPerError - 1];
  //return &states->states[nErrors][states->nStatesPerError[nErrors]-1];
}

void printBacktrace (states_t *states, int depth, int nErrors, size_t stateId) {
  printf("Backtrace:\n");
  state_t *state;
  while ((depth >= 0) || (nErrors > 0)) {
    state = getState(states, depth, nErrors, stateId);
    printf("\tdepth: %i, errors: %i\t", depth, nErrors);
    printState(state, depth); fflush(stdout);
    if (hasTrace(state, MATCH)) {
      --depth;
    }
    else if (hasTrace(state, MISMATCH)) {
      --depth;
      --nErrors;
    }
    else if (hasTrace(state, INSERTION)) {
      --depth;
      assert(nErrors > 0);
      --nErrors;
    }
    else if (hasTrace(state, DELETION)) {
      --nErrors;
    }
    else {
      assert(false);
    }
    stateId = state->previousState;
  }
}

bool isStateInserted (states_t *states, size_t depth, size_t nErrors, bwtinterval_t *interval) {
  //printf("Checking whether interval %" PRIu64 "-%" PRIu64 " is inserted at %zu, %zu, with %zu elements\n", interval->k, interval->l, depth, nErrors, states->nStates[depth][nErrors]); fflush(stdout);
  size_t firstState = states->firstState[depth][nErrors];
  size_t nStates = states->nStates[depth][nErrors];
  int diff;
  if (nStates == 0) {
    return false;
  }
  for (size_t offset = 0; offset < nStates; ) {
    diff = intervalComp(&states->states[nErrors][firstState+offset].interval, interval);
    //printf("  Offset is %zu, state is ", offset);
    //printState(&states->states[nErrors][firstState+offset], depth); fflush(stdout);
    if (diff == 0) {
      //printf("  Got it\n");
      return true;
    }
    if (diff < 0) {
      offset = 2*offset + 1;
    }
    else {
      offset = 2*offset + 2;
    }
  }
  //printf("  Failed\n"); fflush(stdout);
  return false;
}

state_t *insertState (states_t *states, size_t depth, size_t nErrors, bwtinterval_t *interval, unsigned char trace, unsigned char nucleotide, unsigned int previousState, bool simple) {
  //printf("Inserting interval %" PRIu64 "-%" PRIu64 " at %zu, %zu (%c/%c), with %zu elements\n", interval->k, interval->l, depth, nErrors, CIGAR[trace], "ACGT"[nucleotide], states->nStates[depth][nErrors]); fflush(stdout);
  assert((states->nStates[depth][nErrors] == 0) || ((interval->l == 0) == (states->states[nErrors][states->firstState[depth][nErrors]].interval.l == 0)));
  state_t tmpState;
  state_t *state = allocateNewState(states, depth, nErrors);
  state_t *baseState = &states->states[nErrors][states->firstState[depth][nErrors]];
  unsigned int offset = states->nStates[depth][nErrors]-1;
  setState(state, interval, trace, nucleotide, previousState);
  if (simple) return state;
  for (unsigned int parentOffset = (offset-1)/2; (offset != 0) && (stateComp(baseState+parentOffset, baseState+offset) > 0); parentOffset = (offset-1)/2) {
    tmpState = *(baseState+offset); 
    //printf("  Offset is %u (%u), state is ", offset, parentOffset);
    //printState(baseState+offset, depth); fflush(stdout);
    *(baseState+offset) = *(baseState+parentOffset);
    *(baseState+parentOffset) = tmpState;
  } 
  //printf("  Done\n"); fflush(stdout);
  return baseState+offset;
}

state_t *addState (states_t *states, size_t depth, size_t nErrors, bwtinterval_t *interval, unsigned char trace, unsigned char nucleotide, unsigned int previousState, bool add) {
  //state_t *state;
  ++stats->nTentativeStateInsertions;
  //add = true;
  //printf("Trying to add state @ depth %zu, %zu errors, interval %" PRIu64 "-%" PRIu64 " %c/%c -> %u\n", depth, nErrors, interval->k, interval->l, CIGAR[trace], "ACGT"[nucleotide], previousState); fflush(stdout);
  // check if last state has the same interval
  //size_t nStates = states->nStates[depth][nErrors];
  if (! add) {
  //if ((! add) && (nStates > 0)) {
    //nStates += states->firstState[depth][nErrors] - 1;
    //state = &states->states[nErrors][nStates];
    // Taking last state seems more efficient.
    //state = &states->states[nErrors][states->firstState[depth][nErrors]];
    //if ((compareBwtIntervals(&state->interval, interval)) && (state->trace == trace)) {
    if (isStateInserted(states, depth, nErrors, interval)) {
    //if (compareBwtIntervals(&state->interval, interval)) {
      //printf("\t# states: %zu, state id: %zu/%zu-%zu/%zu\n", states->nStates[depth][nErrors], nStates, states->firstState[depth][nErrors], states->nStates[depth][nErrors], states->firstState[depth+1][nErrors]);
      //printStates(states, depth);
      //printf("  Found ");
      //printState(state, depth);
      //fflush(stdout);
      //printBacktrace(states, depth, nErrors, states->nStates[depth][nErrors] - 1);
      ++stats->nSkippedStateInsertions;
      return NULL;
    }
  }
  return insertState(states, depth, nErrors, interval, trace, nucleotide, previousState, add);
  //state = allocateNewState(states, depth, nErrors);
  //setState(state, interval, trace, nucleotide, previousState);
  //if (interval->l != 0) printBacktrace(states, depth, nErrors, states->nStates[depth][nErrors] - 1);
}

states_t *initializeStates(size_t treeSize) {
  states_t *states           = (states_t *)     malloc(sizeof(states_t));
  states->depth              = treeSize + 1 + parameters->maxNErrors;
  states->states             = (state_t **)     malloc((parameters->maxNErrors+1) * sizeof(states_t *));
  states->firstState         = (size_t **)      malloc(states->depth * sizeof(size_t *));
  states->nStates            = (size_t **)      malloc(states->depth * sizeof(size_t *));
  states->allocatedNStates   = (size_t *)       malloc((parameters->maxNErrors+1) * sizeof(size_t));
  states->nStatesPerPosition = (size_t *)       calloc(states->depth,  sizeof(size_t));
  states->minErrors          = (size_t *)       malloc(states->depth * sizeof(size_t));
  states->maxErrors          = (size_t *)       malloc(states->depth * sizeof(size_t));
  states->sw                 = (sw_t *)         malloc(sizeof(sw_t));
  states->bwtBuffer          = (bwt_buffer_t *) malloc(sizeof(bwt_buffer_t));
  for (size_t depth = 0; depth < states->depth; ++depth) {
    states->firstState[depth] = (size_t *) calloc((parameters->maxNErrors+1), sizeof(size_t));
    states->nStates[depth]    = (size_t *) calloc((parameters->maxNErrors+1), sizeof(size_t));
    states->minErrors[depth] = SIZE_MAX;
    states->maxErrors[depth] = SIZE_MAX;
  }
  for (size_t nErrors = 0; nErrors <= parameters->maxNErrors; ++nErrors) {
    states->allocatedNStates[nErrors] = N_STATES;
    states->states[nErrors] = (state_t *) malloc(states->allocatedNStates[nErrors] * sizeof(states_t));
  }
  bwtinterval_t interval = {0, bwt->seq_len};
  addState(states, 0, 0, &interval, 0, 0, 0, true);
  createSW(states->sw, states->depth);
  createBwtBuffer(states->bwtBuffer);
  return states;
}

void backtrackStates(states_t *states, size_t level) {
  //TODO Could be skipped in the levels do not change
  //printf("\t\t\tBacktracking to level %zu\n", level);
  //printStates(states, level);
  for (size_t i = level; i <= states->depth; ++i) {
    if (states->nStatesPerPosition[i] == 0) {
      //printf("\t\t\tFinally2\n");
      //printStates(states, level);
      return;
    }
    states->nStatesPerPosition[i] = 0;
    states->minErrors[i]          = SIZE_MAX;
    states->maxErrors[i]          = SIZE_MAX;
    for (size_t j = 0; j <= parameters->maxNErrors; ++j) {
      states->nStates[i][j] = 0;
      states->firstState[i+1][j] = states->firstState[i][j];
    }
  }
  //printf("\t\t\tFinally\n");
  //printStates(states, level);
}


void computeBacktrace (states_t *states, int depth, int nErrors, size_t stateId, outputSam_t *outputSam) {
  state_t *state;
  char cigar;
  outputSam->backtraceSize = 0;
  while ((depth >= 0) || (nErrors > 0)) {
    state = getState(states, depth, nErrors, stateId);
    //printState(state, 101); fflush(stdout);
    cigar = CIGAR[state->trace];
    if ((outputSam->backtraceSize == 0) || (outputSam->backtraceCigar[outputSam->backtraceSize-1] != cigar)) {
      outputSam->backtraceCigar[outputSam->backtraceSize] = cigar;
      outputSam->backtraceLengths[outputSam->backtraceSize] = 1;
      ++outputSam->backtraceSize;
    }
    else {
      ++outputSam->backtraceLengths[outputSam->backtraceSize-1];
    }
    if (hasTrace(state, MATCH)) {
      --depth;
    }
    else if (hasTrace(state, MISMATCH)) {
      --depth;
      assert(nErrors > 0);
      --nErrors;
    }
    else if (hasTrace(state, INSERTION)) {
      --depth;
      assert(nErrors > 0);
      --nErrors;
    }
    else if (hasTrace(state, DELETION)) {
      --nErrors;
    }
    else {
      assert(false);
    }
    stateId = state->previousState;
  }
}

void freeStates(states_t *states) {
  for (size_t i = 0; i <= parameters->maxNErrors; ++i) {
    free(states->states[i]);
  }
  for (size_t depth = 0; depth < states->depth; ++depth) {
    free(states->firstState[depth]);
    free(states->nStates[depth]);
  }
  free(states->states);
  free(states->nStates);
  free(states->allocatedNStates);
  free(states->firstState);
  free(states->nStatesPerPosition);
  free(states->minErrors);
  free(states->maxErrors);
  free(states->bwtBuffer);
  freeSW(states->sw);
  free(states->sw);
  free(states);
}

#endif
