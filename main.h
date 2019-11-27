#ifndef MAIN_H
#define MAIN_H

#include "constants.h"
#include "helper.h"
#include "stats.h"
#include "parameters.h"
#include "sam.h"
#include "edge.h"
#include "cell.h"
#include "quality.h"
#include "tree.h"
#include "bwt.h"
#include "state.h"
#include "sw.h"
#include "states.h"
#include "shortcut.h"
#include "path.h"


/**
 * Print read in SAM format
 */
void printRead (states_t *states, path_t *path, char *quality, count_t *counts, outputSam_t *outputSam) {
  assert(path->depth <= path->maxDepth);
  size_t depth = path->depth;
  int strand, rid;
  int64_t pos;
  char *seq = path->read + path->readPos;
  char *qual = quality;
  bwtint_t nHits = 0;
  unsigned int nErrors = states->minErrors[depth];
  state_t *theseStates = states->states[depth][nErrors];
  size_t nStates;
  simplifyStates(states, depth, nErrors);
  nStates = states->nStates[depth][nErrors];
  //printf("Quality size: %zu vs %zu\n", strlen(quality), depth);
  //printf("Quality: %s (%p)\n", quality, quality);
  //printf("Read:    %s\n", forwardSeq);
  assert(strlen(quality) == depth);
  //printf("Printing read\n");
  //printPath(path); fflush(stdout);
  //printStates(states, path->depth+1); fflush(stdout);
  writeQname(outputSam, counts);
  outputSam->isBackwardSet = false;
  //printf("depth: %zu/%zu, # states: %zu, nErrors: %u, Seq: %s, qual: %s\n", depth, path->maxDepth, nStates, nErrors, seq, qual); fflush(stdout);
  //printState(&theseStates[nStates-1], path->depth); fflush(stdout);
  memcpy(outputSam->forwardSeq,  seq,  (depth+1) * sizeof(char));
  memcpy(outputSam->forwardQual, qual, (depth+1) * sizeof(char));
  for (size_t i = 0; i < nStates; ++i) {
    nHits += theseStates[i].interval.l - theseStates[i].interval.k + 1;
  }
  if (nHits > parameters->maxNHits) {
    fprintf(outputSam->file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tNH:i:%lu\tNM:i:%u\n", outputSam->qname, seq, qual, nHits, nErrors);
    return;
  }
  if ((nHits == 1) && (nErrors == 0)) {
    pos = bwa_sa2pos(bns, bwt, theseStates[0].interval.k, depth, &strand);
    rid = bns_pos2rid(bns, pos);
    pos = pos - bns->anns[rid].offset + 1;
    printReadUniqueNoError(strand, bns->anns[rid].name, pos, depth, outputSam);
    return;
  }
  //printf("Read: %s, read length: %zu\n", forwardSeq, readLength);
  //printStates(states, depth);
  for (size_t stateId = 0; stateId < nStates; ++stateId) {
    printReadState(&theseStates[stateId], path->depth, nHits, nErrors, outputSam);
  }
}

/**
 * Map without error
 */
bool mapWithoutError (states_t *states, size_t depth, unsigned short nt, size_t nErrors) {
  assert(depth > 0);
  state_t *previousState;
  state_t *nextState;
  bwtinterval_t nextInterval;
  bool mapFound = false;
  //printf("    Mapping %c without error at depth %zu with %zu errors and %zu states\n", "ACGT"[nt], depth, nErrors, states->nStates[depth-1][nErrors]);
  /*
  if (states->nStates[depth-1][nErrors] >= MANY_STATES) {
    states->nStates[depth-1][nErrors] = simplifyStates(states->states[depth-1][nErrors], states->nStates[depth-1][nErrors]);
  }
  */
  for (size_t stateId = 0; stateId < states->nStates[depth-1][nErrors]; ++stateId) {
    previousState = &states->states[depth-1][nErrors][stateId];
    if (goDownBwt(states->bwtBuffer, previousState, nt, &nextInterval)) {
      mapFound = true;
      //printState(&nextState, depth);
      nextState = addState(states, depth, nErrors);
      if (nextState == NULL) {
        //printf("      cannot add state\n");
        return false;
      }
      setState(nextState, &nextInterval, MATCH, previousState);
    }
  }
  //printf("    found map: %s\n", mapFound ? "true" : "false");
  return mapFound;
}

/**
 * Find the mappings with nErrors at depth.
 * Supposes that mappings at depth-1 with nErrors are computed.
 * Supposes that mappings at depth with nErrors-1 are computed.
 */
bool _addError (states_t *states, path_t *path, size_t nErrors, size_t depth) {
  //printf("    depth %zu, %zu errors, nucleotide %c, %zu states, min errors: %zu\n", depth, nErrors, "ACGT"[path->nucleotides[depth-1]],  states->nStates[depth-1][nErrors-1], states->minErrors[depth-1]);
  //printf("      Path is: ");
  //for (size_t i = 0; i < depth; ++i) putchar("ACGT"[path->nucleotides[i]]);
  //putchar('\n');
  //TODO may be skipped?
  bwtinterval_t nextInterval;
  state_t *previousState;
  state_t *nextState;
  if ((states->maxErrors[depth] != SIZE_MAX) && (states->maxErrors[depth] >= nErrors)) {
    //printf("      first case: %zu/%zu/%zu %zu/%i\n", states->maxErrors[depth], nErrors, SIZE_MAX, states->nStates[depth][nErrors], N_STATES);
    return (states->nStates[depth][nErrors] < N_STATES);
  }
  for (size_t stateId = 0; stateId < states->nStates[depth-1][nErrors-1]; ++stateId) {
    previousState = &states->states[depth-1][nErrors-1][stateId];
    //printState(state, path->maxDepth);
    // add insertion
    if (! hasTrace(previousState, DELETION)) {
      nextState = addState(states, depth, nErrors);
      if (nextState == NULL) {
        return false;
      }
      setState(nextState, &previousState->interval, INSERTION, previousState);
    }
    // add mismatches
    for (unsigned short nt = 0; nt < N_NUCLEOTIDES; ++nt) {
      if (goDownBwt(states->bwtBuffer, previousState, nt, &nextInterval)) {
        //addState(states, depth-1, nErrors, &newState);
        if (nt != path->nucleotides[depth-1]) {
          //printState(newState, path->maxDepth);
          nextState = addState(states, depth, nErrors);
          if (nextState == NULL) {
            //printf("      third case\n");
            return false;
          }
          setState(nextState, &nextInterval, MISMATCH | nt, previousState);
        }
      }
    }
  }
  // add deletions
  for (size_t stateId = 0; stateId < states->nStates[depth][nErrors-1]; ++stateId) {
    previousState = &states->states[depth][nErrors-1][stateId];
    if (! hasTrace(previousState, INSERTION)) {
      for (unsigned short nt = 0; nt < N_NUCLEOTIDES; ++nt) {
        if (goDownBwt(states->bwtBuffer, previousState, nt, &nextInterval)) {
          nextState = addState(states, depth, nErrors);
          if (nextState == NULL) {
            return false;
          }
          setState(nextState, &nextInterval, DELETION | nt, previousState);
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
  //TODO adapth this
  /*
  if (states->nStates[depth][nErrors] >= MANY_STATES) {
    states->nStates[depth][nErrors] = simplifyStates(states->states[depth][nErrors], states->nStates[depth][nErrors]);
  }
  */
  //states->nStates[depth][nErrors] = simplifyStates(states->states[depth][nErrors], states->nStates[depth][nErrors]);
  if (states->nStates[depth][nErrors] >= N_STATES) {
    //printf("      fifth case: # states = %zu >= %i after simplification.\n", states->nStates[depth][nErrors], N_STATES);
    return false;
  }
  return true;
}

/**
 * Map the last nucleotide, with one more error, at depth path->depth-1.
 */
bool addError (states_t *states, path_t *path) {
  size_t nErrors = states->minErrors[path->depth-1] + 1;
  size_t firstDepth;
  bool firstDepthFound = false;
  // Find the first place where the mappings with nErrors are computed
  for (firstDepth = path->depth; (firstDepth > 1) && (! firstDepthFound); --firstDepth) {
    for (size_t nucleotide = 0; (nucleotide < N_NUCLEOTIDES) && (! firstDepthFound); ++nucleotide) {
      firstDepthFound = (states->nStates[firstDepth][nErrors] > 0);
      //if (firstDepthFound) printf("Found first depth at depth = %zu, # errors = %zu, nucleotide = %zu\n", firstDepth, nErrors, nucleotide);
    }
  }
  if (firstDepthFound) firstDepth += 2;
  //printf("  adding error with %zu @ %zu from %zu\n", nErrors, path->depth, firstDepth);
  for (size_t depth = firstDepth; depth <= path->depth; ++depth) {
    if (! _addError(states, path, nErrors, depth)) {
      //printf("    too many errors\n");
      return false;
    }
  }
  return true;
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
  //printPath(path);
  //TODO check this
  if (path->depth == TREE_BASE_SIZE) {
    for (size_t nErrors = states->minErrors[TREE_BASE_SIZE]; nErrors <= states->maxErrors[TREE_BASE_SIZE]; ++nErrors) {
      simplifyStates(states, TREE_BASE_SIZE, nErrors);
    }
  }
}

bool shortCutCondition (const states_t *states, const tree_t *tree, const path_t *path) {
  if (path->depth < TREE_BASE_SIZE) return false;
  if (path->edgeLength != 0) return false;
  if (isTerminal(&tree->cells[path->cellIds[path->nCells]])) return false;
  if (! isCellUnbranched(tree, &tree->cells[path->cellIds[path->nCells]])) return false;
  //if (getNChildren(&tree->cells[path->cellIds[path->nCells]]) != 1) return false;
  unsigned int nStates = 0;
  for (unsigned int nErrors = states->minErrors[path->depth]; nErrors <= states->maxErrors[path->depth]; ++nErrors) {
    nStates += states->nStates[path->depth][nErrors];
    if (nStates > MAX_SW_N_STATES) return false;
    for (unsigned int stateId = 0; stateId < states->nStates[path->depth][nErrors]; ++stateId) {
      if (states->states[path->depth][nErrors][stateId].interval.l - states->states[path->depth][nErrors][stateId].interval.k >= MAX_SW_N_SEQUENCES) return false;
    }
  }
  return true;
}

bool tryShortCuts (const tree_t *tree, states_t *states, path_t *path, outputSam_t *outputSam) {
  //printf("    Entering tryShortCuts will cell %" PRIu64 " at depth %zu, last nt %i, and read pos %zu, edge len: %zu\n", path->cellIds[path->nCells], path->depth, path->nucleotides[path->depth-1], path->readPos, path->edgeLength);
  //printPath(path);
  state_t bestStates[2 * MAX_EDGE_LENGTH];
  unsigned int currentNErrors, bestNErrors = parameters->maxNErrors + 1, baseNErrors = 0;
  unsigned int bestStateId = 0;
  unsigned int alignmentSize = 0;
  unsigned int genomeAlignmentSize;
  unsigned int nGenomeSequences;
  unsigned int previousDepth = path->depth;
  unsigned int previousNCells = path->nCells;
  unsigned int previousReadPos = path->readPos;
  uint64_t     bestPos = 0;
  int          bestStrand = 0;
  int          rid;
  char         cigar;
  state_t *previousState;
  state_t *nextState;
  uint64_t cellId = path->cellIds[path->nCells];
  cell_t *cell = &tree->cells[cellId];
  edge_t *edge = getFirstEdge(cell);
  assert(edge != NULL);
  unsetReadSequence(states->sw);
  do {
    assert(edge->length != 0);
    path->edges[path->nCells] = *edge;
    path->cellIds[++path->nCells] = edge->cellId;
    addReadSequence(states->sw, edge->sequence, edge->length);
    cellId = edge->cellId;
    cell = &tree->cells[cellId];
    edge = getFirstEdge(cell);
  }
  while (edge != NULL);
  //TODO store all the best paths instead
  for (unsigned int nErrors = states->minErrors[previousDepth]; (nErrors <= states->maxErrors[previousDepth]) && (nErrors <= bestNErrors); ++nErrors) {
    //printf("      Trying with %u errors\n", nErrors);
    simplifyStates(states, previousDepth, nErrors);
    for (unsigned int stateId = 0; stateId < states->nStates[previousDepth][nErrors]; ++stateId) {
      //printf("      Trying with state #%u/%zu %p\n", stateId, states->nStates[previousDepth][nErrors], &states->nStates[previousDepth][nErrors]);
      //TODO: Do not store all the genome sequences
      nGenomeSequences = setGenomeSequences(states->sw, &states->states[previousDepth][nErrors][stateId]);
      for (unsigned int genomeSequenceId = 0; genomeSequenceId < nGenomeSequences; ++genomeSequenceId) {
        setGenomeSequence(states->sw, genomeSequenceId);
        currentNErrors = getScore(states->sw, genomeSequenceId, bestNErrors - nErrors - 1);
        if (currentNErrors + nErrors < bestNErrors) {
          bestNErrors   = currentNErrors + nErrors;
          baseNErrors   = nErrors;
          bestStateId   = stateId;
          alignmentSize = states->sw->alignmentSize;
          bestPos       = states->sw->poss[genomeSequenceId];
          bestStrand    = states->sw->strands[genomeSequenceId];
          memcpy(bestStates, states->sw->states, alignmentSize * sizeof(state_t));
        }
      }
    }
  }
  if (bestNErrors > parameters->maxNErrors) {
    path->nCells  = previousNCells;
    path->depth   = previousDepth;
    path->readPos = previousReadPos;
    //printf("    Leaving tryShortCuts wrong with cell %" PRIu64 " at depth %zu, last nt %i, and read pos %zu\n", path->cellIds[path->nCells], path->depth, path->nucleotides[path->depth-1], path->readPos);
    return false;
  }
  assert(alignmentSize > 0);
  previousState = &states->states[previousDepth][baseNErrors][bestStateId];
  //printf("\tSize of alignment: %u, starting from %p: %lu, %u, %u\n", alignmentSize, previousState, previousDepth, baseNErrors, bestStateId);
  currentNErrors = baseNErrors;
  genomeAlignmentSize = 0;
  if (bestNErrors == 0) {
    outputSam->backtraceSize       = 1;
    outputSam->backtraceLengths[0] = path->depth;
    outputSam->backtraceCigar[0]   = 'M';
  }
  else {
    computeBacktrace(previousState, outputSam);
  }
  for (unsigned int i = 0; i < alignmentSize; ++i) {
    nextState = addState(states, path->depth, currentNErrors);
    if (nextState == NULL) {
      return false;
    }
    setState(nextState, &bestStates[i].interval, bestStates[i].trace, previousState);
    if ((nextState->trace & BACKTRACE_MASK) != DELETION) {
      //printf("\t\tAdding nucleotide %c\n", "ACGT"[bestStates[i].trace & NUCLEOTIDE_MASK]);
      appendNucleotidePath(path, bestStates[i].trace & NUCLEOTIDE_MASK, DNA5_TO_CHAR[bestStates[i].trace & NUCLEOTIDE_MASK]);
    }
    if (! hasTrace(nextState, MATCH)) {
      ++currentNErrors;
    }
    if (! hasTrace(nextState, INSERTION)) {
      ++genomeAlignmentSize;
    }
    //printf("\tAdding state @ %lu/%u/%lu: ", path->depth, currentNErrors, states->nStates[path->depth][currentNErrors]-1);
    //printState(&states->states[path->depth][currentNErrors][states->nStates[path->depth][currentNErrors]-1], 100);
    previousState = nextState;
    assert(currentNErrors <= bestNErrors);
    cigar = CIGAR[nextState->trace >> BACKTRACE_OFFSET];
    if (bestNErrors != 0) {
      if ((outputSam->backtraceSize == 0) || (outputSam->backtraceCigar[outputSam->backtraceSize-1] != cigar)) {
        outputSam->backtraceCigar[outputSam->backtraceSize] = cigar;
        outputSam->backtraceLengths[outputSam->backtraceSize] = 1;
        ++outputSam->backtraceSize;
      }
      else {
        ++outputSam->backtraceLengths[outputSam->backtraceSize-1];
      }
    }
  }
  assert(currentNErrors == bestNErrors);
  //printStates(states, path->depth+1);
  //printPath(path);
  writeQname(outputSam, cell->counts);
  memcpy(outputSam->forwardSeq,  path->read + path->readPos,  (path->depth+1) * sizeof(char));
  memcpy(outputSam->forwardQual, findQuality(&tree->qualities, cellId), (path->depth+1) * sizeof(char));
  outputSam->isBackwardSet = false;
  computeCigar(outputSam);
  bestPos = (bestStrand)? bestPos + genomeAlignmentSize: bestPos - genomeAlignmentSize;
  rid     = bns_pos2rid(bns, bestPos);
  //printf("Seq: %s, qual: %s, strand: %i, cigar F: %s, cigar B: %s, qname: %s\n", outputSam->forwardSeq, outputSam->forwardQual, bestStrand, outputSam->forwardCigar, outputSam->backwardCigar, outputSam->qname);
  preparePrintRead(bestPos, rid, bestStrand, 1, 1, bestNErrors, outputSam);
  //printf("    Leaving tryShortCuts right with cell %" PRIu64 " at depth %zu, last nt %i, and read pos %zu\n", path->cellIds[path->nCells], path->depth, path->nucleotides[path->depth-1], path->readPos);
  //printPath(path);
  path->nCells  = previousNCells;
  path->depth   = previousDepth;
  path->readPos = previousReadPos;
  backtrackStates(states, path->depth);
  return true;
}

/*
bool tryShortCut (const tree_t *tree, states_t *states, path_t *path, FILE *outputSamFile) {
  return false;
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
*/

void printProgress(path_t *path) {
  if (path->depth == TREE_BASE_SIZE) {
    if ((path->cellIds[path->nCells] & 11111) == 0) {
      fprintf(stderr, "Progress: %" PRIu64 "/%u\n", path->cellIds[path->nCells], N_TREE_BASE);
    }
  }
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
  printProgress(path);
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
void _map (const tree_t *tree, states_t *states, path_t *path, outputSam_t *outputSam) {
  bool mappable = true;
  char *quality;
  while (true) {
    if ((mappable) && (shortCutCondition(states, tree, path))) {
      tryShortCuts(tree, states, path, outputSam);
      mappable = false;
      /*
      if (mappable) {
        //printf("Short cut with positive exit\n"); fflush(stdout);
        if ((quality = findQuality(&tree->qualities, path->cellIds[path->nCells])) != NULL) {
          printRead(states, path, quality, tree->cells[path->cellIds[path->nCells]].counts, outputSam);
        }
      }
      */
    }
    else {
      // advance in the tree
      //printf("Going to next tree\n");
      if (! goNextTree(tree, states, path, mappable)) {
        return;
      }
      //printPath(path); fflush(stdout);
      //printf("%" PRIu64 ": ", path->cellIds[path->nCells]);
      //printCell(&tree->cells[path->cellIds[path->nCells]]);
      //printf("\n");
      // prefix of the read is not mappable
      mappable = findBestMapping(states, path);
      //printf("\tRead is mappable: %s\n", (mappable)? "yes": "no");
      if ((path->depth >= TREE_BASE_SIZE) && (mappable) && (path->edgeLength == 0)) {
        if ((quality = findQuality(&tree->qualities, path->cellIds[path->nCells])) != NULL) {
          printRead(states, path, quality, tree->cells[path->cellIds[path->nCells]].counts, outputSam);
        }
      }
    }
  }
}

/**
 * Allocate/free structures before/after mapping
 */
void map (const tree_t *tree, outputSam_t *outputSam) {
  states_t *states = initializeStates(tree->depth);
  path_t   *path   = initializePath(tree->depth);
  //printf("depth: %zu\n", tree->depth);
  //addState(&states, 0, 0, firstState);
  //goDownTree(tree, &path);
  _map(tree, states, path, outputSam);
  freeStates(states);
  freePath(path);
}

#endif
