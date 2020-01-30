#ifndef SW_H
#define SW_H

#include "constants.h"
#include "helper.h"
#include "state.h"
#include "edge.h"


/******* Smith-Waterman cell type *******/
/**
 * Contains value for backtrace
 */
typedef struct {
  unsigned int value: MAX_SW_COST_SIZE;
  unsigned int backTrace: BACKTRACE_SIZE;
  unsigned int nucleotide: NUCLEOTIDES_BITS;
} sw_cell_t;

/******* Smith-Waterman type *******/
/**
 * Perform the Smith-Waterman algorithm on the reads and the genome.
 * The matrix uses only a band.
 *                  <--- read ---->      original     new
 *                 +-+-+-+-+-+-+-+-+
 *  ^            ^ |X|X|2| | | | | |      +-+-+        +-+
 *  g     # dels | +-+-+-+-+-+-+-+-+      |3|1|        |1|
 *  e            | |X|1|N| | | | | |      +-+-+      +-+-+
 *  n              +-+-+-+-+-+-+-+-+      |2|N|      |3|N|
 *  o   0 error -> |0|N|N| | | | | |      +-+-+      +-+-+
 *  m              +-+-+-+-+-+-+-+-+                 |2|
 *  e            | |1|N|N| | | | | |                 +-+
 *  |      # ins | +-+-+-+-+-+-+-+-+   1: deletion
 *  v            v |2|N|N| | | | | |   2: insertion
 *                 +-+-+-+-+-+-+-+-+   3: (mis)match
 *
 *
 */
typedef struct {
  unsigned short  *readSequence;
  unsigned int     readLength;
  unsigned short **genomeSequences;
  unsigned int    *genomeLengths;
  bwtinterval_t    genomeInterval;
  uint64_t        *poss;
  int             *strands;
  unsigned int     nGenomeSequences;
  unsigned int     alignmentSize;
  sw_cell_t      **matrix;
  state_t         *states;
  size_t           nCols;
} sw_t;

void createSW (sw_t *sw, size_t depth) {
  //printf("SW initialization depth is %zu\n", depth);
  size_t nRows = 2 * parameters->maxNErrors + 1;
  sw->nCols           = depth + 1;
  sw->genomeSequences = (unsigned short **) malloc(MAX_SW_N_SEQUENCES * sizeof(unsigned short *));
  sw->genomeLengths   = (unsigned int *)    calloc(MAX_SW_N_SEQUENCES,  sizeof(unsigned int));
  sw->poss            = (uint64_t *)        malloc(MAX_SW_N_SEQUENCES * sizeof(uint64_t));
  sw->strands         = (int *)             malloc(MAX_SW_N_SEQUENCES * sizeof(int));
  sw->readSequence    = (unsigned short *)  malloc(depth * sizeof(unsigned short));
  for (unsigned int i = 0; i < MAX_SW_N_SEQUENCES; ++i) {
    sw->genomeSequences[i] = (unsigned short *) malloc(depth * sizeof(unsigned short));
  }
  sw->matrix = (sw_cell_t **) malloc(sw->nCols * sizeof(unsigned int *));
  for (size_t i = 0; i < sw->nCols; ++i) {
    sw->matrix[i] = (sw_cell_t *) malloc(nRows * sizeof(unsigned int));
  }
  sw->matrix[0][parameters->maxNErrors].value = 0;
  for (unsigned int i = 1; i <= parameters->maxNErrors; ++i) {
    sw->matrix[0][i+parameters->maxNErrors].value = i;
    sw->matrix[0][i+parameters->maxNErrors].backTrace = SW_DELETION;
  }
  for (unsigned int i = 1; i <= parameters->maxNErrors; ++i) {
    sw->matrix[i][parameters->maxNErrors-i].value = i;
    sw->matrix[i][parameters->maxNErrors-i].backTrace = SW_INSERTION;
  }
  sw->states = (state_t *) malloc(sw->nCols * sizeof(state_t));
  for (unsigned int i = 0; i < sw->nCols; ++i) {
    setEmptyState(&sw->states[i]);
  }
}

void freeSW (sw_t *sw) {
  for (unsigned int i = 0; i < MAX_SW_N_SEQUENCES; ++i) {
    free(sw->genomeSequences[i]);
  }
  for (size_t i = 0; i < sw->nCols; ++i) {
    free(sw->matrix[i]);
  }
  free(sw->genomeSequences);
  free(sw->genomeLengths);
  free(sw->poss);
  free(sw->strands);
  free(sw->readSequence);
  free(sw->matrix);
  free(sw->states);
}

void unsetReadSequence (sw_t *sw) {
  //printf("\t\t\tUnsetting read sequence\n");
  sw->readLength = 0;
}

void addReadSequence (sw_t *sw, sequence_t sequence, unsigned int length) {
  //TODO double-check the direction
  //printf("\t\t\tSetting read sequence (%d): ", length);
  //printSequence(sequence, length);
  //printf(", read length: %d\n", sw->readLength);
  for (unsigned int i = 0; i < length; ++i) {
    sw->readSequence[sw->readLength] = sequence & NUCLEOTIDE_MASK;
    sequence >>= NUCLEOTIDES_BITS;
    if (sw->readLength < parameters->maxNErrors) {
      sw->matrix[sw->readLength+1][parameters->maxNErrors-sw->readLength-1].nucleotide = sw->readSequence[sw->readLength];
    }
    ++sw->readLength;
  }
  //printf("\t\t\tSetting read sequence (%d): ", length);
  //printSequenceLong(sw->readSequence, sw->readLength);
  //printf("\n");
}

unsigned int setGenomeSequences (sw_t *sw, state_t *state) {
  assert(state->interval.k <= state->interval.l);
  int      isRev;
  bwtint_t pos;
  //unsigned int nSequencesUndup = 0;
  sw->nGenomeSequences = state->interval.l - state->interval.k + 1;
  // copy all sequences in the interval
  for (unsigned int sequenceId = 0; sequenceId < sw->nGenomeSequences; ++sequenceId) {
    pos = bwt_sa(bwt, state->interval.k + sequenceId);
    pos = bns_depos(bns, pos, &isRev);
    //printf("Pos: %lu/%lu/%lu\n", pos, bwt->bwt_size, bwt->seq_len); fflush(stdout);
    //bns_pos2rid(bns, pos);
    //printf("\t\t\tSetting genome sequence: ");
    sw->poss[sequenceId]          = pos;
    sw->strands[sequenceId]       = isRev;
    sw->genomeLengths[sequenceId] = 0;
    //for (unsigned int i = 0; i < MAX_SW_LENGTH; ++i) {
    //printf("Max genome length: %zu\n", sw->readLength + parameters->maxNErrors);
    for (unsigned int i = 0; i < sw->readLength + parameters->maxNErrors; ++i) {
      pos = (isRev)? pos+1: pos-1;
      uint64_t c = DNA5_TO_INT_REV[isRev][_get_pac(pac, pos)];
      sw->genomeSequences[sequenceId][i] = c;
      //printf("%c", "ACGT"[c]);
      ++sw->genomeLengths[sequenceId];
      if ((isRev) && (pos == bwt->seq_len)) {
        break;
      }
      else if ((! isRev) && (pos == 0)) {
        break;
      }
    }
    //printf("\n"); fflush(stdout);
    //printf("\t\t\tTransformed into:        ");
    //printSequenceLong(sw->genomeSequences[sequenceId], sw->genomeLengths[sequenceId]);
    //printf("\n"); fflush(stdout);
  }
  sw->genomeInterval = state->interval;
  //TODO remove duplicate sequences
  /*
  for (unsigned int sequenceId = 1; sequenceId < sw->nGenomeSequences; ++sequenceId) {
    if (memcmp(sw->genomeSequences[sequenceId], sw->genomeSequences[nSequencesUndup], sw->genomeLengths[sequenceId] * sizeof(unsigned short)) == 0) {
      sw->genomeSequences[sequenceId] = NULL;
      ++nSequencesUndup;
    }
  }
  */
  return sw->nGenomeSequences;
}

void setGenomeSequence (sw_t *sw, unsigned int genomeSequenceId) {
  //printf("\t\t\tUsing genome sequence: ");
  //printSequenceLong(sw->genomeSequences[genomeSequenceId], sw->genomeLengths[genomeSequenceId]);
  //printf("\n");
  for (unsigned int i = 1; i <= MIN(parameters->maxNErrors, sw->genomeLengths[genomeSequenceId]); ++i) {
    sw->matrix[0][parameters->maxNErrors+i].nucleotide = getNucleotide(sw->genomeSequences[genomeSequenceId][i], i-1);
    //printf("\t\tinsertion (0, %zu): %c\n", parameters->maxNErrors+i, "ACGT"[sw->matrix[0][parameters->maxNErrors+i].nucleotide]);
  }
}

bool tryNoDiffSW (sw_t *sw, unsigned int genomeSequenceId) {
  //printf("\t\t\tShrinking genome to:     ");
  //printSequenceLong(sw->genomeSequences[genomeSequenceId], sw->readLength);
  //printf("\n");
  unsigned short genomeNucleotide;
  if (sw->readLength != sw->genomeLengths[genomeSequenceId]) {
    return false;
  }
  if (memcmp(sw->readSequence, sw->genomeSequences[genomeSequenceId], sw->readLength * sizeof(unsigned short)) != 0) {
    return false;
  }
  sw->alignmentSize = sw->readLength;
  //printf("\t\tFound right away!\n");
  for (unsigned int i = 0; i < sw->genomeLengths[genomeSequenceId]; ++i) {
    genomeNucleotide = sw->genomeSequences[genomeSequenceId][i];
    sw->states[i].trace         = MATCH;
    sw->states[i].nucleotide    = genomeNucleotide;
    sw->states[i].previousState = 0;
    //printf("\t\t\tState @ %u\n", i);
    //printState(&sw->states[i], 100);
    //printf("\n");
  }
  return true;
}

unsigned int getScore (sw_t *sw, unsigned int genomeSequenceId, unsigned int maxNErrors) {
  uint64_t reconstructedGenomeSequence = 0;
  unsigned short readNucleotide, genomeNucleotide;
  unsigned int readId, genomeId, xId, yId, startingYId = 2 * maxNErrors + 2, bestYId = startingYId;
  unsigned int v1, v2, v3;
  unsigned int minValue = maxNErrors + 1;
  unsigned int reversedBacktraceId;
  bwtinterval_t previousInterval;
  bool match;
  if (tryNoDiffSW(sw, genomeSequenceId)) {
    return 0;
  }
  if (maxNErrors == 0) {
    return 1;
  }
  for (unsigned int i = 1; i <= maxNErrors+1; ++i) {
    sw->matrix[0][i].nucleotide = sw->genomeSequences[genomeSequenceId][i-1];
  }
  // local alignment
  //printf("\t\tAt most: %u errors\n", maxNErrors); fflush(stdout);
  for (readId = 0; readId < sw->readLength; ++readId) {
    readNucleotide = sw->readSequence[readId];
    minValue = maxNErrors + 1;
    for (int diff = -maxNErrors; diff <= (int) maxNErrors; ++diff) {
      if (((int) readId) + diff >= 0) {
        genomeId = readId + diff;
        yId = diff + ((int) parameters->maxNErrors);
        if (genomeId < sw->genomeLengths[genomeSequenceId]) {
          genomeNucleotide = sw->genomeSequences[genomeSequenceId][genomeId];
          match = (readNucleotide == genomeNucleotide);
          v1 = (diff == (int) -maxNErrors)? maxNErrors+1: (unsigned int) sw->matrix[readId+1][yId-1].value + 1;
          v2 = (diff == (int) maxNErrors)?  maxNErrors+1: (unsigned int) sw->matrix[readId][yId+1].value + 1;
          v3 = sw->matrix[readId][yId].value + (match? 0: 1);
          //printf("\t\t\t\t\t%u\n", v1); fflush(stdout);
          //printf("\t\t\t\t\t%u\n", v2); fflush(stdout);
          //printf("\t\t\t\t\t%u\n", v3); fflush(stdout);
          //printf("\t\t\t\t(%u, %i/%i) Values %u %u %u / %u\n", readId, diff, yId, v1, v2, v3, maxNErrors); fflush(stdout);
          if ((v3 <= v1) && (v3 <= v2)) {
            sw->matrix[readId+1][yId].value = v3;
            sw->matrix[readId+1][yId].backTrace = (match)? SW_MATCH: SW_MISMATCH;
            sw->matrix[readId+1][yId].nucleotide = genomeNucleotide;
            //printf("\t\t\t\t(mis)match (%u, %u) <- %u (%c)\n", readId+1, yId, v3, "ACGT"[readNucleotide]);
          }
          else if (v1 <= v2) {
            sw->matrix[readId+1][yId].value = v1;
            sw->matrix[readId+1][yId].backTrace = SW_DELETION;
            sw->matrix[readId+1][yId].nucleotide = genomeNucleotide;
            //printf("\t\t\t\tdeletion (%u, %u) <- %u (%c)\n", readId+1, yId, v1, "ACGT"[genomeNucleotide]);
          }
          else {
            sw->matrix[readId+1][yId].value = v2;
            sw->matrix[readId+1][yId].backTrace = SW_INSERTION;
            sw->matrix[readId+1][yId].nucleotide = readNucleotide;
            //printf("\t\t\t\tinsertion (%u, %u) <- %u (%c)\n", readId+1, yId, v2, "ACGT"[readNucleotide]);
          }
          minValue = MIN(minValue, sw->matrix[readId+1][yId].value);
        }
        else {
          sw->matrix[readId+1][yId].value = maxNErrors + 1;
          sw->matrix[readId+1][yId].backTrace = SW_INSERTION;
          sw->matrix[readId+1][yId].nucleotide = readNucleotide;
        }
      }
    }
    //printf("\t\tminValue %u @ %u\n", minValue, readId); fflush(stdout);
    if (minValue == maxNErrors + 1) {
      return maxNErrors + 1;
    }
  }
  //printf("\t\tFound with mismatches\n");
  // get end of backtrace, favoring (mis)matches
  readId = sw->readLength;
  //printf("\t\t\t\tmax #errors %u, readId %u\n", maxNErrors, readId);
  for (unsigned int diff = 0; diff <= maxNErrors; ++diff) {
    //printf("\t\t\t\tdiff %u\n", diff); fflush(stdout);
    if (diff == 0) {
      yId = parameters->maxNErrors;
      //printf("\t\t\t\t\ttest 0: (%u, %u) = %u\n", readId, yId, sw->matrix[readId][yId].value); fflush(stdout);
      assert(sw->matrix[readId][yId].value >= minValue);
      if (sw->matrix[readId][yId].value == minValue) {
        bestYId = yId;
        break;
      }
    }
    else {
      if (readId >= diff) {
        yId = parameters->maxNErrors - diff;
        //printf("\t\t\t\t\ttest -1: (%u, %u) = %u\n", readId, yId, sw->matrix[readId][yId].value); fflush(stdout);
        assert(sw->matrix[readId][yId].value >= minValue);
        if (sw->matrix[readId][yId].value == minValue) {
          bestYId = yId;
          break;
        }
      }
      yId = parameters->maxNErrors + diff;
      //printf("\t\t\t\t\ttest +1: (%u, %u) = %u\n", readId, yId, sw->matrix[readId][yId].value); fflush(stdout);
      assert(sw->matrix[readId][yId].value >= minValue);
      if (sw->matrix[readId][yId].value == minValue) {
        bestYId = yId;
        break;
      }
    }
  }
  //printf("\t\tbest/starting/min: %u/%u/%u\n", bestYId, startingYId, minValue);
  assert(bestYId != startingYId);
  // get size of backtracing
  xId = readId;
  yId = bestYId;
  sw->alignmentSize = 0;
  while((xId != 0) || (yId != parameters->maxNErrors)) {
    ++sw->alignmentSize;
    switch(sw->matrix[xId][yId].backTrace) {
      case SW_MATCH:
      case SW_MISMATCH:
        //printf("\t\t\t\t(mis)match (%i, %i)\n", xId, yId); fflush(stdout);
        assert(xId > 0);
        --xId;
        break;
      case SW_INSERTION:
        //printf("\t\t\t\tinsertion (%i, %i)\n", xId, yId); fflush(stdout);
        assert(xId > 0);
        --xId;
        ++yId;
        break;
      case SW_DELETION:
        //printf("\t\t\t\tdeletion (%i, %i)\n", xId, yId); fflush(stdout);
        assert(yId > 0);
        --yId;
        break;
      default:
        assert(false);
    }
  }
  //printf("\t\t\tdone, size: %u\n", sw->alignmentSize); fflush(stdout);
  // backtrace
  xId = readId;
  yId = bestYId;
  for (unsigned int backtraceId = 0; backtraceId < sw->alignmentSize; ++backtraceId) {
    reversedBacktraceId = sw->alignmentSize - backtraceId - 1;
    sw->states[reversedBacktraceId].trace      = sw->matrix[xId][yId].backTrace;
    sw->states[reversedBacktraceId].nucleotide = sw->matrix[xId][yId].nucleotide;
    switch(sw->matrix[xId][yId].backTrace) {
      case SW_MATCH:
      case SW_MISMATCH:
        //printf("\t\t\t\t(mis)match (%u, %u) %c\n", xId, yId, "ACGT"[sw->matrix[xId][yId].nucleotide]); fflush(stdout);
        reconstructedGenomeSequence <<= NUCLEOTIDES_BITS;
        reconstructedGenomeSequence |= sw->matrix[xId][yId].nucleotide;
        --xId;
        break;
      case SW_INSERTION:
        //printf("\t\t\t\tinsertion (%u, %u) %c\n", xId, yId, "ACGT"[sw->matrix[xId][yId].nucleotide]); fflush(stdout);
        --xId;
        ++yId;
        break;
      case SW_DELETION:
        //printf("\t\t\t\tdeletion (%u, %u)\n", xId, yId); fflush(stdout);
        reconstructedGenomeSequence <<= NUCLEOTIDES_BITS;
        reconstructedGenomeSequence |= sw->matrix[xId][yId].nucleotide;
        --yId;
        break;
      default:
        assert(false);
    }
  }
  assert(xId == 0);
  assert(yId == parameters->maxNErrors);
  // add BWT intervals
  previousInterval = sw->genomeInterval;
  for (unsigned int backtraceId = 0; backtraceId < sw->alignmentSize; ++backtraceId) {
  //for (int backtraceId = sw->alignmentSize - 1; backtraceId >= 0; --backtraceId) {
    //printf("\t\t\tbacktrace: %u\n", sw->states[backtraceId].trace);
    if (! hasTrace(&sw->states[backtraceId], INSERTION)) {
      genomeNucleotide = reconstructedGenomeSequence & NUCLEOTIDE_MASK;
      /*
      bwt_2occ(bwt, previousInterval.k-1, previousInterval.l, genomeNucleotide, &sw->states[backtraceId].interval.k, &sw->states[backtraceId].interval.l);
      sw->states[backtraceId].interval.k = bwt->L2[genomeNucleotide] + sw->states[backtraceId].interval.k + 1;
      sw->states[backtraceId].interval.l = bwt->L2[genomeNucleotide] + sw->states[backtraceId].interval.l;
      */
      sw->states[backtraceId].previousState = 0;
      reconstructedGenomeSequence >>= NUCLEOTIDES_BITS;
    }
    else {
      sw->states[backtraceId].interval = previousInterval;
      sw->states[backtraceId].previousState = 0;
    }
    //printf("\t\t\tState @ %u\n", backtraceId);
    //printState(&sw->states[backtraceId], 100);
    //printf("\n"); fflush(stdout);
    previousInterval = sw->states[backtraceId].interval;
  }
  //printf("\t\t\tgenome sequence: %lu\n", reconstructedGenomeSequence); fflush(stdout);
  assert(reconstructedGenomeSequence == 0);
  return minValue;
}

#endif
