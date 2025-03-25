#ifndef STATES_H
#define STATES_H

#include "constants.h"
#include "parameters.h"
#include "state.h"
#include "bwt.h"

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
  state_t      **states;
  size_t       **firstState;
  size_t       **nStates;
  size_t        *allocatedNStates;
  size_t        *nStatesPerPosition;
  size_t        *minErrors;
  size_t        *maxErrors;
  size_t        depth;
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

/**
 * Update counts and allocate arrays to receive new states
 */
void updateCounts (states_t *states, size_t depth, size_t nErrors, size_t nStates) {
  //printf("Entering update counts with states %p depth %zu, %zu errors, and %zu states, errors [%zu-%zu/%lu]\n", states, depth, nErrors, nStates, states->minErrors[depth], states->maxErrors[depth], SIZE_MAX); fflush(stdout);
  assert(depth <= states->depth);
  size_t nStatesPerError;
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
  }
  /*
  if (depth > 0) {
    printf("Previous state: %zu/%zu\n", states->firstState[depth-1][nErrors], states->nStates[depth-1][nErrors]);
    fflush(stdout);
  }
  printf("first state: %zu\n", states->firstState[depth][nErrors]); fflush(stdout);
  */
  states->nStates[depth][nErrors]   += nStates;
  states->nStatesPerPosition[depth] += nStates;
  nStatesPerError = states->firstState[depth][nErrors] + states->nStates[depth][nErrors];
  //stats->maxNStates[nErrors] = MAX(stats->maxNStates[nErrors], nStatesPerError);
  if (nStatesPerError >= states->allocatedNStates[nErrors]-1) {
    while (nStatesPerError >= states->allocatedNStates[nErrors]-1) {
      states->allocatedNStates[nErrors] *= 2;
    }
    states->states[nErrors] = (state_t *) reallocOrDie(states->states[nErrors], states->allocatedNStates[nErrors] * sizeof(state_t));
    //printf("Reallocating states with %zu errors to %zu states\n", nErrors, states->allocatedNStates[nErrors]);
    if (states->states[nErrors] == NULL) {
      fprintf(stderr, "Memory error: Cannot allocate an array of States for errors %zu of size %zu * %zu\n", nErrors, states->allocatedNStates[nErrors], sizeof(state_t));
      exit(EXIT_FAILURE);
    }
    //printf("Reallocation #states with %zu errors: now %zu\n", nErrors, states->allocatedNStates[nErrors]); fflush(stdout);
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

state_t *_addState (states_t *states, size_t depth, size_t nErrors, bwtinterval_t *interval, unsigned char trace, unsigned char nucleotide, unsigned int previousState) {
  //printf("Inserting interval %" PRIu64 "-%" PRIu64 " at %zu, %zu (%c/%c), with %zu elements\n", interval->k, interval->l, depth, nErrors, CIGAR[trace], "ACGT"[nucleotide], states->nStates[depth][nErrors]); fflush(stdout);
  assert((states->nStates[depth][nErrors] == 0) || ((interval->l == 0) == (states->states[nErrors][states->firstState[depth][nErrors]].interval.l == 0)));
  state_t *state;
  updateCounts(states, depth, nErrors, 1);
  state = &states->states[nErrors][states->firstState[depth][nErrors] + states->nStates[depth][nErrors] - 1];
  setState(state, interval, trace, nucleotide, previousState);
  //printf("Inserted interval %" PRIu64 "-%" PRIu64 " at %zu, %zu (%c/%c), with %zu elements\n", interval->k, interval->l, depth, nErrors, CIGAR[trace], "ACGT"[nucleotide], states->nStates[depth][nErrors]); fflush(stdout);
  //printStates(states, depth); fflush(stdout);
  return state;
}

size_t getNIndels (states_t *states, int depth, int nErrors, size_t stateId) {
  state_t *state;
  size_t nIndels = 0;
  while ((depth > 0) || (nErrors > 0)) {
    state = getState(states, depth, nErrors, stateId);
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
      ++nIndels;
    }
    else if (hasTrace(state, DELETION)) {
      assert(nErrors > 0);
      --nErrors;
      ++nIndels;
    }
    else {
      assert(false);
    }
    stateId = state->previousState;
  }
  return nIndels;
}

void putBestStateFirst (states_t *states, size_t depth, size_t nErrors, size_t firstStateId, size_t lastStateId) {
  state_t *theseStates = states->states[nErrors] + states->firstState[depth][nErrors];
  size_t minNIndels = getNIndels(states, depth, nErrors, firstStateId);
  size_t minId      = firstStateId;
  for (size_t stateId = firstStateId + 1; (stateId < lastStateId) && (minNIndels > 0); ++stateId) {
    size_t nIndels = getNIndels(states, depth, nErrors, stateId);
    if (nIndels < minNIndels) {
      minNIndels = nIndels;
      minId = stateId;
    }
  }
  if (minId != firstStateId) {
    theseStates[firstStateId] = theseStates[minId];
  }
}

// Aim: discard redundant mappings (favors mismatches over indels)
void simplifyStates (states_t *states, size_t depth, size_t nErrors) {
  assert(depth < states->depth);
  assert(nErrors <= states->maxErrors[depth]);
  size_t nStates = states->nStates[depth][nErrors];
  state_t *theseStates = states->states[nErrors] + states->firstState[depth][nErrors];
  size_t firstStateId = 0;
  size_t secondStateId = 1;
  size_t thisNIndels = getNIndels(states, depth, nErrors, 0);
  size_t thatNIndels = 0;
  if (nStates <= 1) {
    return;
  }
  //TODO check that the pointers are good!
  //printf("Entering Simplify States @ depth %zu with %zu errors and %zu elements.\n", depth, nErrors, nStates); fflush(stdout);
  // Sort the states by position
  qsort(theseStates, nStates, sizeof(state_t), sortCompareStates);
  //printStates(states, depth);
  /* ==BEFORE==
   *   STEP 1         |   ...   |              STEP 2        |   ...   |       STEP 3        |   ...   |
   *                  |   ...   |   e s                      |   ...   |                     |   ...   |
   * nextNState ----> | state x |-\ q t     nextNState ----> | state x |                     | state x |
   *                  |   ...   | | u a                      | state y |    nextNState ----> | state y |
   *                  |   ...   |-/ i t                      |   ...   |                     |   ...   |
   * secondStateId -> | state y |   v e     secondStateId -> | state y |                     | state y |
   *                  |   ...   |   . s                      |   ...   |    secondStateId -> |   ...   |
   * 
   * ==NOW==
   *   STEP 1         |   ...   |              STEP 2        |   ...   |       STEP 3        |   ...   |
   *                  |   ...   |   e s                      |   ...   |                     |   ...   |
   * firstStateId --> | state 1 |-\ q t     firstStateId --> | state 1 |                     | state 1 |
   *                  | state 2 | | u a                      | state 4 |    firstStateId --> | state 4 |
   *                  | state 3 |-/ i t                      | state 5 |    secondStateId -> | state 5 |
   * secondStateId -> | state 4 |   v e     secondStateId -> |   ...   |                     |   ...   |
   *                  | state 5 |   . s
   *                  |   ...   |
   */
  //printf("Starting with depth %zu and %zu errors\n", depth, nErrors);
  //printf("Before:\n"); for (size_t tmp = 0; tmp < nStates; ++tmp) printState(&theseStates[tmp], depth);
  //printStates(states, depth);
  for (; secondStateId < nStates; ++secondStateId) {
    //printf("\tCurrent state: %zu/%zu/%zu\n", firstStateId, secondStateId, nStates); fflush(stdout);
    if (areStatesEqual(&theseStates[firstStateId], &theseStates[secondStateId])) { 
      thatNIndels = getNIndels(states, depth, nErrors, secondStateId);
      if (thatNIndels < thisNIndels) {
        theseStates[firstStateId] = theseStates[secondStateId];
        thisNIndels = thatNIndels;
      }
    }
    else {
      ++firstStateId;
      theseStates[firstStateId] = theseStates[secondStateId];
      thisNIndels = getNIndels(states, depth, nErrors, firstStateId);
      secondStateId = firstStateId; // Will be increased at the end of the for loop
      assert(firstStateId <= nStates);
      assert(secondStateId < N_STATES);
    }
    //printf("\tCurrent state (2): %zu/%zu/%zu\n", firstStateId, secondStateId, nStates); fflush(stdout);
  }
  nStates = firstStateId + 1;
  //printf("to %zu\n", nStates);
  //printf("After\n"); for (size_t tmp = 0; tmp < nStates; ++tmp) printState(&theseStates[tmp], depth);
  //printStates(states, depth);
  assert(states->nStatesPerPosition[depth] >= nStates);
  states->nStatesPerPosition[depth] = states->nStatesPerPosition[depth] - states->nStates[depth][nErrors] + nStates;
  states->nStates[depth][nErrors] = nStates;
}

state_t *addState (states_t *states, size_t depth, size_t nErrors, bwtinterval_t *interval, unsigned char trace, unsigned char nucleotide, unsigned int previousState) {
  return _addState(states, depth, nErrors, interval, trace, nucleotide, previousState);
}

void clearStates (states_t *states) {
  for (size_t depth = 0; depth < states->depth; ++depth) {
    memset(states->firstState[depth], 0, (parameters->maxNErrors+1) * sizeof(size_t));
    memset(states->nStates[depth],    0, (parameters->maxNErrors+1) * sizeof(size_t));
    states->minErrors[depth] = SIZE_MAX;
    states->maxErrors[depth] = SIZE_MAX;
  }
  bwtinterval_t interval = {0, bwt->seq_len};
  addState(states, 0, 0, &interval, 0, 0, 0);
}

states_t *initializeStates(size_t treeSize) {
  states_t *states           = (states_t *)     mallocOrDie(sizeof(states_t));
  states->depth              = treeSize + 1 + parameters->maxNErrors;
  states->states             = (state_t **)     mallocOrDie((parameters->maxNErrors+1) * sizeof(states_t *));
  states->firstState         = (size_t **)      mallocOrDie((states->depth+1) * sizeof(size_t *)); // Allocate one unit more for the backtrack.
  states->nStates            = (size_t **)      mallocOrDie(states->depth * sizeof(size_t *));
  states->allocatedNStates   = (size_t *)       mallocOrDie((parameters->maxNErrors+1) * sizeof(size_t));
  states->nStatesPerPosition = (size_t *)       callocOrDie(states->depth+1,  sizeof(size_t)); // Allocate one unit more for the backtrack.
  states->minErrors          = (size_t *)       mallocOrDie(states->depth * sizeof(size_t));
  states->maxErrors          = (size_t *)       mallocOrDie(states->depth * sizeof(size_t));
  for (size_t depth = 0; depth < states->depth; ++depth) {
    states->firstState[depth] = (size_t *) mallocOrDie((parameters->maxNErrors+1) * sizeof(size_t));
    states->nStates[depth]    = (size_t *) mallocOrDie((parameters->maxNErrors+1) * sizeof(size_t));
  }
  states->firstState[states->depth] = (size_t *) callocOrDie(parameters->maxNErrors+1, sizeof(size_t));
  states->allocatedNStates[0] = treeSize + 3;
  states->states[0] = (state_t *) malloc(states->allocatedNStates[0] * sizeof(state_t));
  if (states->states[0] == NULL) {
    fprintf(stderr, "Error! Cannot allocate sufficient memory (%zu cells) for states with 0 errors.\nExiting.\n", states->allocatedNStates[0]);
    exit(EXIT_FAILURE);
  }
  for (size_t nErrors = 1; nErrors <= parameters->maxNErrors; ++nErrors) {
    states->allocatedNStates[nErrors] = N_STATES;
    states->states[nErrors] = (state_t *) malloc(states->allocatedNStates[nErrors] * sizeof(state_t));
    if (states->states[nErrors] == NULL) {
      fprintf(stderr, "Error! Cannot allocate sufficient memory (%zu cells) for states with %zu errors.\nExiting.\n", states->allocatedNStates[nErrors], nErrors);
      exit(EXIT_FAILURE);
    }
  }
  clearStates(states);
  //createBwtBuffer(&states->bwtBuffer);
  return states;
}

void backtrackStates(states_t *states, size_t level) {
  //printf("\t\t\tBacktracking to level %zu\n", level);
  //printStates(states, level);
  for (size_t i = level; i <= states->depth; ++i) {
    //printf("\t\t\t\tnow %zu/%zu -> %zu\n", i, states->depth, states->nStatesPerPosition[i]); fflush(stdout);
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


size_t computeBacktrace (states_t *states, int depth, int nErrors, size_t stateId, outputSam_t *outputSam) {
  state_t *state;
  char cigar;
  int readSize = 0;
  outputSam->backtraceSize = 0;
  outputSam->nIndels = 0;
  //printf("Compute BT\n");
  //while ((depth >= 0) || (nErrors > 0)) {
  while ((depth > 0) || (nErrors > 0)) {
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
      ++readSize;
    }
    else if (hasTrace(state, MISMATCH)) {
      --depth;
      assert(nErrors > 0);
      --nErrors;
      ++readSize;
    }
    else if (hasTrace(state, INSERTION)) {
      --depth;
      assert(nErrors > 0);
      --nErrors;
      ++outputSam->nIndels;
    }
    else if (hasTrace(state, DELETION)) {
      --nErrors;
      ++readSize;
      ++outputSam->nIndels;
    }
    else {
      assert(false);
    }
    stateId = state->previousState;
  }
  return readSize;
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
  free(states);
}

#endif
