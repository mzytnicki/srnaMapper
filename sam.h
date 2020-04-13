#ifndef SAM_H
#define SAM_H

#include "constants.h"
#include "helper.h"
#include "parameters.h"

#define INITIAL_SAM_SIZE 0x10000
#define MAX_SAM_LINE     0x200

typedef struct {
  FILE *file;
  char *qname;
  unsigned int backtraceSize;
  char *backtraceCigar;
  unsigned char *backtraceLengths;
  char *forwardCigar;
  char *backwardCigar;
  char *forwardSeq;
  char *backwardSeq;
  char *forwardQual;
  char *backwardQual;
  bool isBackwardSet;
  pthread_mutex_t *writeMutex;
  char *output;
  size_t sizeOutput;
  size_t allocatedOutput;
} outputSam_t;

void createOutputSam (outputSam_t *outputSam, unsigned int size, pthread_mutex_t *mutex) {
  outputSam->isBackwardSet    = false;
  outputSam->qname            = (char *)          malloc(MAX_QNAME_SIZE * sizeof(char));
  outputSam->backtraceCigar   = (char *)          malloc(MAX_CIGAR_SIZE * sizeof(char));
  outputSam->backtraceLengths = (unsigned char *) malloc(MAX_CIGAR_SIZE * sizeof(unsigned char));
  outputSam->forwardCigar     = (char *)          malloc(MAX_CIGAR_SIZE * sizeof(char));
  outputSam->backwardCigar    = (char *)          malloc(MAX_CIGAR_SIZE * sizeof(char));
  outputSam->forwardSeq       = (char *)          malloc((size+1)       * sizeof(char));
  outputSam->backwardSeq      = (char *)          malloc((size+1)       * sizeof(char));
  outputSam->forwardQual      = (char *)          malloc((size+1)       * sizeof(char));
  outputSam->backwardQual     = (char *)          malloc((size+1)       * sizeof(char));
  outputSam->writeMutex       = mutex;
  outputSam->sizeOutput       = 0;
  outputSam->allocatedOutput  = INITIAL_SAM_SIZE;
  outputSam->output           = (char *)          malloc(outputSam->allocatedOutput * sizeof(char));
  outputSam->output[0]        = 0;
}

void computeCigar (outputSam_t *outputSam) {
  size_t forwardCigarSize = 0;
  size_t backwardCigarSize = 0;
  //TODO Optimize this
  forwardCigarSize = backwardCigarSize = 0;
  for (size_t backtraceId = 0; backtraceId < outputSam->backtraceSize; ++backtraceId) {
    forwardCigarSize  += sprintf(outputSam->forwardCigar+forwardCigarSize,   "%i%c", outputSam->backtraceLengths[backtraceId], outputSam->backtraceCigar[backtraceId]);
    backwardCigarSize += sprintf(outputSam->backwardCigar+backwardCigarSize, "%i%c", outputSam->backtraceLengths[outputSam->backtraceSize-backtraceId-1], outputSam->backtraceCigar[outputSam->backtraceSize-backtraceId-1]);
  }
  outputSam->forwardCigar[forwardCigarSize] = 0;
  outputSam->backwardCigar[backwardCigarSize] = 0;
}

void setCigarNoError (outputSam_t *outputSam, unsigned int readLength) {
  sprintf(outputSam->forwardCigar, "%uM", readLength);
  sprintf(outputSam->backwardCigar, "%uM", readLength);
}

void computeReverseComplement (outputSam_t *outputSam) {
  size_t readLength = strlen(outputSam->forwardSeq);
  outputSam->isBackwardSet = true;
  reverseComplementSequence(outputSam->forwardSeq, outputSam->backwardSeq, readLength);
  reverseSequence(outputSam->forwardQual, outputSam->backwardQual, readLength);
}

void freeOutputSam (outputSam_t *outputSam) {
  free(outputSam->qname);
  free(outputSam->backtraceCigar);
  free(outputSam->backtraceLengths);
  free(outputSam->forwardCigar);
  free(outputSam->backwardCigar);
  free(outputSam->forwardSeq);
  free(outputSam->backwardSeq);
  free(outputSam->forwardQual);
  free(outputSam->backwardQual);
  free(outputSam->output);
}

void writeQname (outputSam_t *outputSam, count_t *counts) {
  size_t qnameLength = sprintf(outputSam->qname, "read%lu_x", ++nReads);
  for (size_t readsFileId = 0; readsFileId < parameters->nReadsFiles; ++readsFileId) {
    qnameLength += sprintf(outputSam->qname+qnameLength, "_%" PRIu32, counts[readsFileId]);
  }
}

unsigned int computeMapq (unsigned int nHits, unsigned int nErrors) {
  if ((nHits > 1) || (nErrors >= 40)) {
    return 0;
  }
  else {
    return 40 - nErrors;
  }
}

void updateOutputAllocation (outputSam_t *outputSam) {
  if (outputSam->sizeOutput + MAX_SAM_LINE >= outputSam->allocatedOutput) {
    outputSam->allocatedOutput *= 2;
    outputSam->output = (char *) realloc(outputSam->output, outputSam->allocatedOutput * sizeof(char)); 
  }
}

/**
 * Print one line in the SAM file, too many hits
 */
void printReadLineManyHits (char *seq, char *qual, bwtint_t nHits, unsigned int nErrors, outputSam_t *outputSam) {
  outputSam->sizeOutput += sprintf(outputSam->output + outputSam->sizeOutput, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tNH:i:%lu\tNM:i:%u\n", outputSam->qname, seq, qual, nHits, nErrors);
  updateOutputAllocation(outputSam);
}

/**
 * Print one line in the SAM file
 */
void printReadLine (bool forward, unsigned int flag, char *chrName, int64_t pos, bwtint_t nHits, bwtint_t hitId, unsigned int nErrors, outputSam_t *outputSam) {
  char *cigar, *seq, *qual;
  if (! outputSam->isBackwardSet) {
    computeReverseComplement(outputSam);
  }
  cigar = (forward)? outputSam->forwardCigar: outputSam->backwardCigar;
  seq   = (forward)? outputSam->forwardSeq:   outputSam->backwardSeq;
  qual  = (forward)? outputSam->forwardQual:  outputSam->backwardQual;
  outputSam->sizeOutput += sprintf(outputSam->output + outputSam->sizeOutput, "%s\t%u\t%s\t%" PRId64 "\t%d\t%s\t*\t0\t0\t%s\t%s\tNH:i:%lu\tHI:i:%lu\tIH:i:%lu\tNM:i:%u\n", outputSam->qname, flag, chrName, pos, computeMapq(nHits, nErrors), cigar, seq, qual, nHits, hitId, nHits, nErrors);
  updateOutputAllocation(outputSam);
}

/**
 * Print a read, which maps only once, and with no error
 */
void printReadUniqueNoError (int strand, char *chrName, int64_t pos, size_t readLength, outputSam_t *outputSam) {
  unsigned int flag = 0;
  setCigarNoError(outputSam, readLength);
  if (strand == 0) {
    flag = CIGAR_REVERSE;
  }
  printReadLine(strand != 0, flag, chrName, pos, 1, 1, 0, outputSam);
  //fprintf(outputSamFile, "%s\t%u\t%s\t%" PRId64 "\t40\t%zuM\t*\t0\t0\t%s\t%s\tNH:i:1\tHI:i:1\tIH:i:1\tNM:i:0\n", qname, flag, chrName, pos, readLength, seq, qual);
}

/**
 * Collect information before printing read
 */
void preparePrintRead (uint64_t pos, int rid, int strand, bwtint_t nHits, bwtint_t hitId, unsigned int nErrors, outputSam_t *outputSam) {
  unsigned int flag = (hitId == 0)? 0: CIGAR_SECONDARY_HIT;
  if (strand == 0) {
    flag |= CIGAR_REVERSE;
  }
  printReadLine(strand != 0, flag, bns->anns[rid].name, pos, nHits, hitId, nErrors, outputSam);
}

void writeToSam (outputSam_t *outputSam) {
  //printf("Acquiring w mutex\n"); fflush(stdout);
  if (pthread_mutex_lock(outputSam->writeMutex) != 0) {
    fprintf(stderr, "Error, Cannot acquire write mutex!\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  //printf("Acquiring w mutex done\n"); fflush(stdout);
  fprintf(outputSam->file, "%s", outputSam->output);
  //printf("Releasing w mutex\n"); fflush(stdout);
  if (pthread_mutex_unlock(outputSam->writeMutex) != 0) {
    fprintf(stderr, "Error, Cannot release write mutex!\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  //printf("Releasing w mutex done\n"); fflush(stdout);
  outputSam->output[0] = 0;
  outputSam->sizeOutput = 0;
}

FILE *openSamFile() {
  FILE *outputSamFile = fopen(parameters->outputSamFileName, "w");
  if (outputSamFile == NULL) {
    return NULL;
  }
  fprintf(outputSamFile, "@HD\tVN:1.6\tSO:unsorted\n");
  for (int i = 0; i < bns->n_seqs; ++i) {
    fprintf(outputSamFile, "@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
  }
  return outputSamFile;
}


#endif
