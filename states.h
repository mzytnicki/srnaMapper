#ifndef STATES_H
#define STATES_H

#include "constants.h"
#include "parameters.h"
#include "state.h"
#include "bwt.h"
#include "sw.h"

typedef struct {
  size_t nElements;
  size_t size;
  char  *bitField;
} bit_array_t;

void initializeBitArray (bit_array_t *bitArray, size_t size) {
  bitArray->nElements = ((size + (CHAR_BIT) - 1) / (CHAR_BIT));
  bitArray->bitField = (char *) calloc(bitArray->nElements, sizeof(char));
  bitArray->size = size;
}

void freeBitArray (bit_array_t *bitArray) {
  free(bitArray->bitField);
}


void setBitArray (bit_array_t *bitArray, size_t value) {
  assert(value < bitArray->size);
  bitArray->bitField[value / CHAR_BIT] |= (1 << (value % CHAR_BIT));
}

bool isSetBitArray(const bit_array_t *bitArray, size_t value) {
  assert(value < bitArray->size);
  return ((bitArray->bitField[value / CHAR_BIT] & (1 << (value % CHAR_BIT)))? true: false);
}

void resetBitArray (bit_array_t *bitArray) {
  memset(bitArray->bitField, 0, bitArray->nElements * sizeof(char));
}


typedef struct {
  state_t     *statesHash;
  size_t      *idUsed;
  bit_array_t  isUsed;
  size_t       nHashUsed;
  state_t     *statesVector;
  size_t       nStatesVectorAllocated;
  size_t       nStatesVectorUsed;
} states_hash_t;

void initializeStatesHash (states_hash_t *statesHash) {
  statesHash->statesHash             = (state_t *) malloc(N_STATES_HASH_SIZE * sizeof(state_t));
  statesHash->idUsed                 = (size_t *)   malloc(N_STATES_HASH_SIZE * sizeof(size_t));
  statesHash->nHashUsed              = 0;
  statesHash->statesVector           = (state_t *) malloc(N_STATES_VECTOR_SIZE * sizeof(state_t));
  statesHash->nStatesVectorAllocated = N_STATES_VECTOR_SIZE;
  statesHash->nStatesVectorUsed      = 0;
  initializeBitArray(&statesHash->isUsed, N_STATES_HASH_SIZE);
}

void freeStatesHash (states_hash_t *statesHash) {
  free(statesHash->statesHash);
  free(statesHash->idUsed);
  free(statesHash->statesVector);
  freeBitArray(&statesHash->isUsed);
}

state_t *addStateToHashVector (states_hash_t *statesHash, bwtinterval_t *interval, unsigned char trace, unsigned char nucleotide, unsigned int previousState) {
  state_t *state;
  for (size_t stateId = 0; stateId < statesHash->nStatesVectorUsed; ++stateId) {
    if (compareBwtIntervals(interval, &statesHash->statesVector[stateId].interval)) {
      ++stats->nSkippedStateInsertions;
      return &statesHash->statesVector[stateId];
    }
  }
  if (statesHash->nStatesVectorUsed == statesHash->nStatesVectorAllocated-1) {
    statesHash->nStatesVectorAllocated *= 2;
    statesHash->statesVector = (state_t *) realloc(statesHash->statesVector, statesHash->nStatesVectorAllocated);
  }
  state = &statesHash->statesVector[statesHash->nStatesVectorUsed];
  setState(state, interval, trace, nucleotide, previousState);
  ++statesHash->nStatesVectorUsed;
  stats->maxHashVectorSize = MAX(stats->maxHashVectorSize, statesHash->nStatesVectorUsed);
  return state;
}

state_t *_addStateToHash (states_hash_t *statesHash, size_t hashValue, bwtinterval_t *interval, unsigned char trace, unsigned char nucleotide, unsigned int previousState) {
  state_t *state = &statesHash->statesHash[hashValue];
  setBitArray(&statesHash->isUsed, hashValue);
  setState(state, interval, trace, nucleotide, previousState);
  statesHash->idUsed[statesHash->nHashUsed] = hashValue;
  ++statesHash->nHashUsed;
  //stats->maxHashSize = MAX(stats->maxHashSize, statesHash->nHashUsed);
  return state;
}

state_t *addStateToHash (states_hash_t *statesHash, bwtinterval_t *interval, unsigned char trace, unsigned char nucleotide, unsigned int previousState) {
  size_t hashValue = getBwtHash(*interval) % N_STATES_HASH_SIZE;
  if (isSetBitArray(&statesHash->isUsed, hashValue)) {
    if (compareBwtIntervals(interval, &statesHash->statesHash[hashValue].interval)) {
      ++stats->nSkippedStateInsertions;
      return &statesHash->statesHash[hashValue];
    }
    return addStateToHashVector(statesHash, interval, trace, nucleotide, previousState);
  }
  return _addStateToHash(statesHash, hashValue, interval, trace, nucleotide, previousState);
}

void clearHash (states_hash_t *statesHash) {
  statesHash->nHashUsed         = 0;
  statesHash->nStatesVectorUsed = 0;
  resetBitArray(&statesHash->isUsed);
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
  state_t      **states;
  size_t       **firstState;
  size_t       **nStates;
  size_t        *allocatedNStates;
  size_t        *nStatesPerPosition;
  size_t        *minErrors;
  size_t        *maxErrors;
  size_t        depth;
  sw_t          *sw;
  bwt_buffer_t   bwtBuffer;
  states_hash_t  statesHash;
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
  //printf("Entering update counts with depth %zu, %zu errors, and %zu states\n", depth, nErrors, nStates); fflush(stdout);
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
  stats->maxNStates[nErrors] = MAX(stats->maxNStates[nErrors], nStatesPerError);
  if (nStatesPerError >= states->allocatedNStates[nErrors]-1) {
    while (nStatesPerError >= states->allocatedNStates[nErrors]-1) {
      states->allocatedNStates[nErrors] *= 2;
    }
    states->states[nErrors] = (state_t *) realloc(states->states[nErrors], states->allocatedNStates[nErrors] * sizeof(state_t));
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

void addHashStates (states_t *states, size_t depth, size_t nErrors) {
  assert(depth <= states->depth);
  size_t nNewStates = states->statesHash.nHashUsed + states->statesHash.nStatesVectorUsed;
  size_t offset;
  if (nNewStates == 0) {
    return;
  }
  updateCounts(states, depth, nErrors, nNewStates);
  offset = states->firstState[depth][nErrors];
  //printf("Adding %zu+%zu=%zu hash states with %zu errors, starting point: %zu/%p\n", states->statesHash.nHashUsed, states->statesHash.nStatesVectorUsed, nNewStates, nErrors, offset, &states->states[nErrors][offset]); fflush(stdout);
  for (size_t stateId = 0; stateId < states->statesHash.nHashUsed; ++stateId, ++offset) {
    //printf("  Setting hash interval #%zu/%zu %" PRIu64 "-%" PRIu64 "\n", stateId, states->statesHash.idUsed[stateId], states->statesHash.statesHash[states->statesHash.idUsed[stateId]].interval.k, states->statesHash.statesHash[states->statesHash.idUsed[stateId]].interval.l);
    memcpy(&states->states[nErrors][offset], &states->statesHash.statesHash[states->statesHash.idUsed[stateId]], sizeof(state_t));
  }
  for (size_t stateId = 0; stateId < states->statesHash.nStatesVectorUsed; ++stateId, ++offset) {
    memcpy(&states->states[nErrors][offset], &states->statesHash.statesVector[stateId], sizeof(state_t));
  }
  clearHash(&states->statesHash);
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
  return state;
}

state_t *addState (states_t *states, size_t depth, size_t nErrors, bwtinterval_t *interval, unsigned char trace, unsigned char nucleotide, unsigned int previousState, bool directAdd) {
  ++stats->nTentativeStateInsertions;
  if (directAdd) {
    return _addState(states, depth, nErrors, interval, trace, nucleotide, previousState);
  }
  return addStateToHash(&states->statesHash, interval, trace, nucleotide, previousState);
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
  for (size_t depth = 0; depth < states->depth; ++depth) {
    states->firstState[depth] = (size_t *) calloc((parameters->maxNErrors+1), sizeof(size_t));
    states->nStates[depth]    = (size_t *) calloc((parameters->maxNErrors+1), sizeof(size_t));
    states->minErrors[depth] = SIZE_MAX;
    states->maxErrors[depth] = SIZE_MAX;
  }
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
  bwtinterval_t interval = {0, bwt->seq_len};
  addState(states, 0, 0, &interval, 0, 0, 0, true);
  createSW(states->sw, states->depth);
  createBwtBuffer(&states->bwtBuffer);
  initializeStatesHash(&states->statesHash);
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
  freeStatesHash(&states->statesHash);
  freeSW(states->sw);
  free(states->sw);
  free(states);
}

#endif
