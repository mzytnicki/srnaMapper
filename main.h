#ifndef MAIN_H
#define MAIN_H

#include "constants.h"
#include "helper.h"
#include "parameters.h"
#include "sam.h"
#include "edge.h"
#include "cell.h"
#include "quality.h"
#include "tree.h"
#include "bwt.h"
#include "state.h"
#include "states.h"
#include "path.h"


/**
 * Print read in SAM format
 */
void printRead (states_t *states, path_t *path, count_t *counts, char *quality, char *readNames, outputSam_t *outputSam) {
  assert(path->depth <= path->maxDepth);
  assert(quality != NULL);
  assert(counts != NULL);
  //printf("Entering 'printRead'\n");
  //printCellInfo(cellInfo);
  size_t depth = path->depth;
  int strand, rid;
  int64_t pos;
  char *seq = path->read + path->readPos;
  char *qual = quality;
  bwtint_t hitId = 0;
  bwtint_t nHits = 0;
  unsigned int nErrors = states->minErrors[depth];
  state_t *theseStates = getState(states, depth, nErrors, 0);
  size_t nStates;
  simplifyStates(states, depth, nErrors);
  nStates = states->nStates[depth][nErrors];
  //printf("Quality size: %zu vs %zu\n", strlen(quality), depth);
  //printf("Quality: %s (%p)\n", quality, quality);
  //printf("Read:    %s\n", seq);
  assert(strlen(quality) == depth);
  //printf("Printing read\n");
  //printPath(path); fflush(stdout);
  //printStates(states, path->depth+1); fflush(stdout);
  setCountsSam(outputSam, counts);
  outputSam->isBackwardSet = false;
  //printf("depth: %zu/%zu, # states: %zu, nErrors: %u, Seq: %s, qual: %s\n", depth, path->maxDepth, nStates, nErrors, seq, qual); fflush(stdout);
  //printState(&theseStates[nStates-1], path->depth); fflush(stdout);
  //printf("depth: %zu/%zu, # states: %zu, nErrors: %u, Seq: %s, qual: %s\n", depth, path->maxDepth, nStates, nErrors, seq, qual); fflush(stdout);
  //printState(&theseStates[nStates-1], path->depth); fflush(stdout);
  /*
  for (size_t i = 0; i < path->depth; ++i) {
    printState(getState(states, i, 0, 0), i);
  }
  */
  outputSam->readNames = readNames;
  memcpy(outputSam->forwardSeq,  seq,  (depth+1) * sizeof(char));
  memcpy(outputSam->forwardQual, qual, (depth+1) * sizeof(char));
  for (size_t i = 0; i < nStates; ++i) {
    nHits += theseStates[i].interval.l - theseStates[i].interval.k + 1;
  }
  if (nHits > parameters->maxNHits) {
    printReadLineManyHits(seq, qual, nHits, nErrors, outputSam);
    return;
  }
  if ((nHits == 1) && (nErrors == 0)) {
    pos = bwa_sa2pos(bns, bwt, theseStates[0].interval.k, depth, &strand);
    rid = bns_pos2rid(bns, pos);
    pos = pos - bns->anns[rid].offset + 1;
    printReadUniqueNoError(strand, rid, pos, depth, outputSam);
    updateReadIds(outputSam);
    return;
  }
  //printf("Read: %s, read length: %zu\n", forwardSeq, readLength);
  //printStates(states, depth);
  for (size_t stateId = 0; stateId < nStates; ++stateId) {
    if (nErrors == 0) {
      setCigarNoError(outputSam, path->depth);
    }
    else {
      //printf("depth: %zu\n", depth);
      computeBacktrace(states, path->depth, nErrors, stateId, outputSam);
      computeCigar(outputSam);
    }
    printReadState(&theseStates[stateId], path->depth, nHits, hitId, nErrors, outputSam);
    hitId += theseStates[stateId].interval.l - theseStates[stateId].interval.k + 1;
  }
  removeDuplicatesOutputLines(outputSam, nHits);
  updateReadIds(outputSam);
}

/**
 * Map without error
 */
bool mapWithoutError (states_t *states, size_t depth, unsigned short nt, size_t nErrors /*, bool oneInsertion */) {
  assert(depth > 0);
  state_t *previousState;
  //state_t *nextState;
  bwtinterval_t nextInterval;
  bool mapFound = false;
  //printf("    Mapping %c without error at depth %zu with %zu errors and %zu states (%p) \n", "ACGT"[nt], depth, nErrors, states->nStates[depth-1][nErrors], states); fflush(stdout);
  //printStates(states, depth); fflush(stdout);
  /*
  if (states->nStates[depth-1][nErrors] >= MANY_STATES) {
    simplifyStates(states, depth-1, nErrors);
  }
  */
  for (size_t stateId = 0; stateId < states->nStates[depth-1][nErrors]; ++stateId) {
    previousState = getState(states, depth-1, nErrors, stateId);
    //printf("Previous state with %zu errors, @ depth %zu, id %zu: %p\n", nErrors, depth-1, stateId, previousState);
    //printState(previousState, depth);
    //printStates(states, depth); fflush(stdout);
    //++stats->nBwtPerDepth[depth-1];
    if (goDownBwt(/* &states->bwtBuffer, */ previousState, nt, &nextInterval)) {
      mapFound = true;
      //nextState = addState(states, depth, nErrors);
      //setState(nextState, &nextInterval, MATCH, 0, stateId);
      addState(states, depth, nErrors, &nextInterval, MATCH, nt, stateId);
      //printf("Next state is %p\n", nextState);
      //printState(nextState, depth);
      //printStates(states, depth+1); fflush(stdout);
    }
  }
  //printf("    found map: %s\n", mapFound ? "true" : "false"); fflush(stdout);
  /*
  if (oneInsertion) {
    addHashStates(states, depth, nErrors);
  }
  */
  return mapFound;
}

/**
 * Find the mappings with nErrors at depth.
 * Return true if at least a state was added.
 * Supposes that mappings at depth-1 with nErrors are computed.
 * Supposes that mappings at depth with nErrors-1 are computed.
 */
bool _addError (states_t *states, path_t *path, size_t nErrors, size_t depth) {
  //printf("    depth %zu, %zu errors, nucleotide %c, %zu states, min errors: %zu\n", depth, nErrors, "ACGT"[path->nucleotides[depth-1]],  states->nStates[depth-1][nErrors-1], states->minErrors[depth-1]);
  //printf("      Path is: ");
  //for (size_t i = 0; i < depth; ++i) putchar("ACGT"[path->nucleotides[i]]);
  //putchar('\n');
  //TODO may be skipped?
  //bwtinterval_t nextIntervals[N_NUCLEOTIDES];
  bwtinterval_t nextInterval;
  state_t *previousState;
  bool stateAdded = false;
  //state_t *nextState;
  if ((states->maxErrors[depth] != SIZE_MAX) && (states->maxErrors[depth] >= nErrors)) {
    //printf("      first case: %zu/%zu/%zu %zu/%i\n", states->maxErrors[depth], nErrors, SIZE_MAX, states->nStates[depth][nErrors], N_STATES);
    return false;
  }
  for (size_t stateId = 0; stateId < states->nStates[depth-1][nErrors-1]; ++stateId) {
    previousState = getState(states, depth-1, nErrors-1, stateId);
    //printState(state, path->maxDepth);
    // add insertion
    if (! hasTrace(previousState, DELETION)) {
      addState(states, depth, nErrors, &previousState->interval, INSERTION, 0, stateId);
      stateAdded = true;
      //nextState = addState(states, depth, nErrors);
      //setState(nextState, &previousState->interval, INSERTION, 0, stateId);
    }
    // add mismatches
    /*
    goDownBwt3Nt(states->bwtBuffer, previousState, path->nucleotides[depth-1], nextIntervals);
    for (unsigned short nt = 0; nt < N_NUCLEOTIDES; ++nt) {
      ++stats->nBwtPerDepth[depth-1];
      if (nt != path->nucleotides[depth-1]) {
    */
    for (unsigned short nt = 0; nt < N_NUCLEOTIDES; ++nt) {
      //++stats->nBwtPerDepth[depth-1];
      if (nt != path->nucleotides[depth-1]) {
        //addState(states, depth-1, nErrors, &newState);
        if (goDownBwt(/* &states->bwtBuffer, */ previousState, nt, &nextInterval)) {
          //printState(newState, path->maxDepth);
          //nextState = addState(states, depth, nErrors);
          //setState(nextState, &nextInterval, MISMATCH, nt, stateId);
          addState(states, depth, nErrors, &nextInterval, MISMATCH, nt, stateId);
          stateAdded = true;
        }
      }
    }
  }
  // add deletions
  for (size_t stateId = 0; stateId < states->nStates[depth][nErrors-1]; ++stateId) {
    previousState = getState(states, depth, nErrors-1, stateId);
    if (! hasTrace(previousState, INSERTION)) {
      //goDownBwt4Nt(states->bwtBuffer, previousState, nextIntervals);
      for (unsigned short nt = 0; nt < N_NUCLEOTIDES; ++nt) {
        //++stats->nBwtPerDepth[depth];
        if (goDownBwt(/* &states->bwtBuffer, */ previousState, nt, &nextInterval)) {
          //nextState = addState(states, depth, nErrors);
          //setState(nextState, &nextInterval, DELETION, nt, stateId);
          addState(states, depth, nErrors, &nextInterval, DELETION, nt, stateId);
          stateAdded = true;
        }
      }
    }
  }
  // add match
  if (mapWithoutError(states, depth, path->nucleotides[depth-1], nErrors)) {
    stateAdded = true;
  }
  //TODO adapth this
  /*
  if (states->nStates[depth][nErrors] >= MANY_STATES) {
    states->nStates[depth][nErrors] = simplifyStates(states->states[depth][nErrors], states->nStates[depth][nErrors]);
  }
  */
  //states->nStates[depth][nErrors] = simplifyStates(states->states[depth][nErrors], states->nStates[depth][nErrors]);
  //addHashStates(states, depth, nErrors);
  return stateAdded;
}

/**
 * Map the last nucleotide, with one more error, at depth path->depth-1.
 */
bool addError (states_t *states, path_t *path) {
  size_t nErrors = states->minErrors[path->depth-1] + 1;
  size_t firstDepth;
  bool firstDepthFound = false;
  bool stateAdded = false;
  // Find the first place where the mappings with nErrors are computed
  for (firstDepth = path->depth; (firstDepth > 1) && (! firstDepthFound); --firstDepth) {
    //for (size_t nucleotide = 0; (nucleotide < N_NUCLEOTIDES) && (! firstDepthFound); ++nucleotide) {
      firstDepthFound = (states->nStates[firstDepth][nErrors] > 0);
      //if (firstDepthFound) printf("Found first depth at depth = %zu, # errors = %zu, nucleotide = %zu\n", firstDepth, nErrors, nucleotide);
    //}
  }
  if (firstDepthFound) firstDepth += 2;
  //printf("  adding error with %zu @ %zu from %zu\n", nErrors, path->depth, firstDepth);
  for (size_t depth = firstDepth; depth <= path->depth; ++depth) {
    if (_addError(states, path, nErrors, depth)) {
      //printf("    too many errors\n");
      stateAdded = true;
    }
  }
  return stateAdded;
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
  //printPath(path); fflush(stdout);
  //TODO check this
  if (path->depth == TREE_BASE_SIZE) {
    for (size_t nErrors = MAX(1, states->minErrors[TREE_BASE_SIZE]); nErrors <= states->maxErrors[TREE_BASE_SIZE]; ++nErrors) {
      simplifyStates(states, TREE_BASE_SIZE, nErrors);
    }
  }
  else if (states->nStates[path->depth][states->minErrors[path->depth]] >= MANY_STATES) {
    simplifyStates(states, path->depth, states->minErrors[path->depth]);
  }
}

bool shortCutCondition (const states_t *states, const tree2_t *tree, const path_t *path) {
  state_t *state;
  if (path->depth < TREE_BASE_SIZE) return false;
  if (path->edgeLength != 0) return false;
  if (isTerminal(&tree->cells[path->cellIds[path->nCells]])) return false;
  if (! isCell2Unbranched(tree, path->cellIds[path->nCells])) return false;
  //if (getNChildren(&tree->cells[path->cellIds[path->nCells]]) != 1) return false;
  unsigned int nStates = 0;
  for (unsigned int nErrors = states->minErrors[path->depth]; nErrors <= states->maxErrors[path->depth]; ++nErrors) {
    nStates += states->nStates[path->depth][nErrors];
    if (nStates > MAX_SW_N_STATES) return false;
    for (unsigned int stateId = 0; stateId < states->nStates[path->depth][nErrors]; ++stateId) {
      state = getState(states, path->depth, nErrors, stateId);
      if (state->interval.l - state->interval.k >= MAX_SW_N_SEQUENCES) return false;
    }
  }
  return true;
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
  //printProgress(path);
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
void _map (const tree2_t *tree, states_t *states, path_t *path, cellVisitor_t *cellVisitor, uint32_t firstCellId, uint32_t lastCellId, outputSam_t *outputSam) {
  bool mappable = true;
  while (true) {
    //if (path->depth > TREE_BASE_SIZE) printPath(path);
    //if (path->depth > TREE_BASE_SIZE) printStates(states, path->depth);
    //printf("%" PRIu32 ": ", path->cellIds[path->nCells]);
    //printCell2(&tree->cells[path->cellIds[path->nCells]]);
    //printf("\n");
    //if (path->depth > TREE_BASE_SIZE) printf("Going to next tree\n");
    if (! goNextTree2(tree, states, path, firstCellId, lastCellId, mappable)) {
      return;
    }
    // prefix of the read is not mappable
    mappable = findBestMapping(states, path);
    //if (path->depth > TREE_BASE_SIZE) printf("\tRead is mappable: %s\n", (mappable)? "yes": "no");
    if ((path->depth >= TREE_BASE_SIZE) && (mappable) && (path->edgeLength == 0)) {
      //if (path->depth > TREE_BASE_SIZE) printPath(path);
      //if ((cellInfo = getCellInfoTree(path->cellIds[path->nCells], cellVisitor)) != NULL) {
      if (getCellInfoTree(tree, path->cellIds[path->nCells], cellVisitor)) {
        printRead(states, path, getCounts(&tree->cellInfos, *cellVisitor), getQuality(&tree->cellInfos, *cellVisitor), getCellInfoReadNames(&tree->cellInfos, *cellVisitor), outputSam);
      }
    }
  }
  writeToSam(outputSam, true);
}

#endif
