#ifndef STATE_H
#define STATE_H

#include <limits.h>
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
  unsigned char trace:      BACKTRACE_SIZE;
  unsigned char nucleotide: NUCLEOTIDES_BITS;
  unsigned int previousState;
};

void setEmptyState (state_t *state) {
  state->interval.k = state->interval.l = 0;
  state->trace         = 0;
  state->nucleotide    = 0;
  state->previousState = -1;
}

void setState (state_t *state, bwtinterval_t *interval, unsigned char trace, unsigned char nucleotide, unsigned int previousState) {
  //printf("Setting state to %" PRIu64 "-%" PRIu64 ", %c\n", interval->k, interval->l, "ACGT"[nucleotide]);
  state->interval      = *interval;
  state->trace         = trace;
  state->nucleotide    = nucleotide;
  state->previousState = previousState;
}

bwtinterval_t *getStateInterval (state_t *state) {
  return &state->interval;
}

bool hasTrace (state_t *state, unsigned char traceType) {
  return (state->trace == traceType);
}

void printState (state_t *state, size_t maxDepth) {
  //assert(getStateInterval(state)->k <= getStateInterval(state)->l);
  size_t s = 0;
  char *tmpSeq1, *tmpSeq2;
  tmpSeq1 = (char *) mallocOrDie(255);
  tmpSeq2 = (char *) mallocOrDie(255);
  tmpSeq1[maxDepth] = 0;
  tmpSeq2[maxDepth] = 0;
  int is_rev;
  if (state->interval.l == 0) {
    bwtint_t p1 = bwt_sa(bwt, getStateInterval(state)->k); // position on the forward-reverse coordinate
    p1 = bns_depos(bns, p1, &is_rev); // position on the forward strand; this may be the first base or the last base
    for (size_t i = 0; i < maxDepth; ++i) {
      bwtint_t p = is_rev? p1-i: p1+i;
      if (p == ULONG_MAX) {
        maxDepth = i;
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
        maxDepth = i;
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
  }
  tmpSeq2[s] = 0;
  int      strand;
  int64_t  pos = bwa_sa2pos(bns, bwt, getStateInterval(state)->k, maxDepth, &strand);
  int      rid = bns_pos2rid(bns, pos);
  pos          = pos - bns->anns[rid].offset + 1;
  printf("\t\t\t\t%" PRId64 "-%" PRIu64 " (%p, %s:%" PRIu64 " (%i)): %s (%c/%c -> %u)\n", getStateInterval(state)->k, getStateInterval(state)->l, state, bns->anns[rid].name, pos, strand, tmpSeq2, CIGAR[state->trace], "ACGT"[state->nucleotide], state->previousState);
}

bool areStatesEqual (state_t *state1, state_t *state2) {
  return ((getStateInterval(state1)->k == getStateInterval(state2)->k) && (getStateInterval(state1)->l == getStateInterval(state2)->l));
}

int sortCompareStates (const void *state1, const void *state2) {
  return (getStateInterval((state_t *) state1)->k - (getStateInterval((state_t *) state2)->k));
}

/**
 * Print all the reads in a given state
 */
void printReadState (state_t *state, size_t depth, bwtint_t nHits, unsigned int hitId, unsigned int nErrors, outputSam_t *outputSam) {
  int64_t pos;
  int strand, rid;
  for (bwtint_t bwtint = getStateInterval(state)->k; bwtint <= getStateInterval(state)->l; ++bwtint) {
    ++hitId;
    pos = bwa_sa2pos(bns, bwt, bwtint, depth, &strand);
    rid = bns_pos2rid(bns, pos);
    pos = pos - bns->anns[rid].offset + 1;
    preparePrintRead(pos, rid, strand, nHits, hitId, nErrors, outputSam);
  }
}


#endif
