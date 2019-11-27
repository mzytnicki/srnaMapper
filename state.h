#ifndef STATE_H
#define STATE_H

#include "constants.h"
#include "sam.h"

/******* State type *******/
/**
 * A state is an prefix of the genome tree, with additional information to perform backtrace.
 * It is:
 *   - a BWT-interval
 *   - the current operation of the backtrace
 *   - the parent state for a backtrace
 */

typedef struct state_t state_t;

struct state_t {
  bwtinterval_t interval;
  unsigned char trace;
  struct state_t *previousState;
};

void setState (state_t *state, bwtinterval_t *interval, unsigned char trace, state_t *previousState) {
  state->interval      = *interval;
  state->trace         = trace;
  state->previousState = previousState;
}

bwtinterval_t *getStateInterval (state_t *state) {
  return &state->interval;
}

bool hasTrace (state_t *state, unsigned char traceType) {
  return ((state->trace & BACKTRACE_MASK) == traceType);
}

void printState (state_t *state, size_t maxDepth) {
  //assert(getStateInterval(state)->k <= getStateInterval(state)->l);
  size_t s = 0;
  char *tmpSeq1, *tmpSeq2;
  tmpSeq1 = (char *) malloc(255);
  tmpSeq2 = (char *) malloc(255);
  tmpSeq1[maxDepth] = 0;
  tmpSeq2[maxDepth] = 0;
  int is_rev;
  bwtint_t p1 = bwt_sa(bwt, getStateInterval(state)->k); // position on the forward-reverse coordinate
  p1 = bns_depos(bns, p1, &is_rev); // position on the forward strand; this may be the first base or the last base
  for (size_t i = 0; i < maxDepth; ++i) {
    bwtint_t p = is_rev? p1-i: p1+i;
    if (p == ULONG_MAX) {
      tmpSeq1[i] = 0;
      break;
    }
    tmpSeq1[i] = ((is_rev)? "TGCAN": "ACGTN")[_get_pac(pac, p)];
  }
  bwtint_t p2 = bwt_sa(bwt, getStateInterval(state)->l); // position on the forward-reverse coordinate
  p2 = bns_depos(bns, p2, &is_rev); // position on the forward strand; this may be the first base or the last base
  for (size_t i = 0; i < maxDepth; ++i) {
    bwtint_t p = is_rev? p2-i: p2+i;
    if (p == ULONG_MAX) {
      tmpSeq2[i] = 0;
      break;
    }
    tmpSeq2[i] = ((is_rev)? "TGCAN": "ACGTN")[_get_pac(pac, p)];
  }
  for (s = 0; s < maxDepth; ++s) {
    if (tmpSeq1[s] != tmpSeq2[s]) {
      tmpSeq1[s] = 0;
      break;
    }
  }
  for (size_t i = 0; i < s; ++i) {
    tmpSeq2[i] = tmpSeq1[s-i-1];
  }
  tmpSeq2[s] = 0;
  printf("\t\t\t\t%" PRIu64 "-%" PRIu64 " (%p): %s (%c/%c -> %p)\n", getStateInterval(state)->k, getStateInterval(state)->l, state, tmpSeq2, CIGAR[(state->trace & BACKTRACE_MASK) >> BACKTRACE_OFFSET], "ACGT"[state->trace & NUCLEOTIDE_MASK], state->previousState);
}

void computeBacktrace (state_t *state, outputSam_t *outputSam) {
  assert(state != NULL);
  outputSam->backtraceSize = 0;
  while (state->previousState != NULL) {
    //printState(state, 101); fflush(stdout);
    char cigar = CIGAR[state->trace >> BACKTRACE_OFFSET];
    if ((outputSam->backtraceSize == 0) || (outputSam->backtraceCigar[outputSam->backtraceSize-1] != cigar)) {
      outputSam->backtraceCigar[outputSam->backtraceSize] = cigar;
      outputSam->backtraceLengths[outputSam->backtraceSize] = 1;
      ++outputSam->backtraceSize;
    }
    else {
      ++outputSam->backtraceLengths[outputSam->backtraceSize-1];
    }
    state = state->previousState;
  }
}

/**
 * Print all the reads in a given state
 */
void printReadState (state_t *state, size_t depth, bwtint_t nHits, unsigned int nErrors, outputSam_t *outputSam) {
  unsigned int hitId = 0;
  int64_t pos;
  int strand, rid;
  if (nErrors != 0) {
    //printf("depth: %zu\n", depth);
    computeBacktrace(state, outputSam);
    computeCigar(outputSam);
  }
  for (bwtint_t bwtint = getStateInterval(state)->k; bwtint <= getStateInterval(state)->l; ++bwtint) {
    ++hitId;
    pos = bwa_sa2pos(bns, bwt, bwtint, depth, &strand);
    rid = bns_pos2rid(bns, pos);
    pos = pos - bns->anns[rid].offset + 1;
    preparePrintRead(pos, rid, strand, nHits, hitId, nErrors, outputSam);
  }
}


#endif
