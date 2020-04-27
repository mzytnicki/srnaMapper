#ifndef SAM_H
#define SAM_H

#include "constants.h"
#include "helper.h"
#include "parameters.h"

#define INITIAL_SAM_SIZE 0x10000
#define MAX_SAM_LINE     0x200

#define MAX_CHR_SIZE 255

/**
 * Stores a hit, and the number of corresponding identical sequences in the input files
 */
typedef struct {
  count_t     *counts;
  unsigned int flag;
  char         chr[MAX_CHR_SIZE];
  uint64_t     pos;
  char         cigar[MAX_CIGAR_SIZE];
  char        *sequence;
  char        *quality;
  unsigned int nHits;
  unsigned int hitId;
  unsigned int mapq;
  unsigned int nErrors;
} samLine_t;

void createSamLine (samLine_t *samLine, size_t maxReadLength, unsigned int nCounts) {
  samLine->counts   = (count_t *) malloc(nCounts * sizeof(count_t));
  samLine->sequence = (char *)    malloc(maxReadLength+1);
  samLine->quality  = (char *)    malloc(maxReadLength+1);
}

void freeSamLine (samLine_t *samLine) {
  free(samLine->counts);
  free(samLine->sequence);
  free(samLine->quality);
}

unsigned int computeMapq (unsigned int nHits, unsigned int nErrors) {
  if ((nHits > 1) || (nErrors >= 40)) {
    return 0;
  }
  else {
    return 40 - nErrors;
  }
}

/**
 * Copy one SAM information to the buffer
 */
void copyToSamLine (samLine_t *samLine, size_t nCounts, count_t *counts, unsigned int flag, char *chrName, int64_t pos, bwtint_t nHits, bwtint_t hitId, unsigned int nErrors, char *cigar, char *sequence, char *quality) {
  memcpy(samLine->counts,   counts,   nCounts * sizeof(count_t));
  strcpy(samLine->chr,      chrName);
  strcpy(samLine->cigar,    cigar);
  strcpy(samLine->sequence, sequence);
  strcpy(samLine->quality,  quality);
  samLine->flag    = flag;
  samLine->pos     = pos;
  samLine->nHits   = nHits;
  samLine->hitId   = hitId;
  samLine->nErrors = nErrors;
  samLine->mapq    = computeMapq(nHits, nErrors);
} 


void printSamTailLine (FILE *file, samLine_t *samLine) {
  fprintf(file, "\t%u\t%s\t%" PRId64 "\t%d\t%s\t*\t0\t0\t%s\t%s\tNH:i:%u\tHI:i:%u\tIH:i:%u\tNM:i:%u\n", samLine->flag, samLine->chr, samLine->pos, samLine->mapq, samLine->cigar, samLine->sequence, samLine->quality, samLine->nHits, samLine->hitId, samLine->nHits, samLine->nErrors);
}

void printSamLineUnique (FILE *file, samLine_t *samLine) {
  fprintf(file, "read%lu", ++nReads);
  for (size_t readsFileId = 0; readsFileId < parameters->nReadsFiles; ++readsFileId) {
    fprintf(file, "_%" PRIu32, samLine->counts[readsFileId]);
  }
  printSamTailLine(file, samLine);
}

void printSamLineMultiple (FILE *file, samLine_t *samLine, unsigned int fileId, unsigned long int readId) {
  for (unsigned int countId = 0; countId < samLine->counts[fileId]; ++countId) {
    fprintf(file, "read%lu", readId);
    printSamTailLine(file, samLine);
  }
}

#define N_SAM_LINES 0x100

typedef struct {
  FILE            **outputFiles;
  count_t          *counts;
  unsigned int      sequenceSize;
  unsigned int      backtraceSize;
  char             *backtraceCigar;
  unsigned char    *backtraceLengths;
  char             *forwardCigar;
  char             *backwardCigar;
  char             *forwardSeq;
  char             *backwardSeq;
  char             *forwardQual;
  char             *backwardQual;
  bool             isBackwardSet;
  pthread_mutex_t *writeMutex;
  samLine_t       *samLines;
  size_t           nSamLines;
} outputSam_t;

void createOutputSam (outputSam_t *outputSam, unsigned int readSize, pthread_mutex_t *mutex) {
  outputSam->sequenceSize     = readSize;
  outputSam->isBackwardSet    = false;
  outputSam->counts           = (count_t *)       malloc(parameters->nReadsFiles * sizeof(count_t));
  outputSam->backtraceCigar   = (char *)          malloc(MAX_CIGAR_SIZE          * sizeof(char));
  outputSam->backtraceLengths = (unsigned char *) malloc(MAX_CIGAR_SIZE          * sizeof(unsigned char));
  outputSam->forwardCigar     = (char *)          malloc(MAX_CIGAR_SIZE          * sizeof(char));
  outputSam->backwardCigar    = (char *)          malloc(MAX_CIGAR_SIZE          * sizeof(char));
  outputSam->forwardSeq       = (char *)          malloc((readSize+1)            * sizeof(char));
  outputSam->backwardSeq      = (char *)          malloc((readSize+1)            * sizeof(char));
  outputSam->forwardQual      = (char *)          malloc((readSize+1)            * sizeof(char));
  outputSam->backwardQual     = (char *)          malloc((readSize+1)            * sizeof(char));
  outputSam->writeMutex       = mutex;
  outputSam->nSamLines        = 0;
  outputSam->samLines         = (samLine_t *) malloc(N_SAM_LINES * sizeof(samLine_t));
  for (size_t samLineId = 0; samLineId < N_SAM_LINES; ++samLineId) {
    createSamLine(&outputSam->samLines[samLineId], readSize, parameters->nReadsFiles);
  }
}

void freeOutputSam (outputSam_t *outputSam) {
  free(outputSam->counts);
  free(outputSam->backtraceCigar);
  free(outputSam->backtraceLengths);
  free(outputSam->forwardCigar);
  free(outputSam->backwardCigar);
  free(outputSam->forwardSeq);
  free(outputSam->backwardSeq);
  free(outputSam->forwardQual);
  free(outputSam->backwardQual);
  for (size_t samLineId = 0; samLineId < N_SAM_LINES; ++samLineId) {
    freeSamLine(&outputSam->samLines[samLineId]);
  }
  free(outputSam->samLines);
}

/**
 * Create CIGAR string from the backtrace
 */
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

/**
 * Compute the reverse-complement of the sequence and the quality
 */
void computeReverseComplement (outputSam_t *outputSam) {
  size_t readLength = strlen(outputSam->forwardSeq);
  outputSam->isBackwardSet = true;
  reverseComplementSequence(outputSam->forwardSeq, outputSam->backwardSeq, readLength);
  reverseSequence(outputSam->forwardQual, outputSam->backwardQual, readLength);
}

/**
 * Print all the SAM buffers to the output file
 */
void writeToSam (outputSam_t *outputSam) {
  if (pthread_mutex_lock(outputSam->writeMutex) != 0) {
    fprintf(stderr, "Error, Cannot acquire write mutex!\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  if (parameters->uniqueOutputFile) {
    for (unsigned int samLineId = 0; samLineId < outputSam->nSamLines; ++samLineId) {
      printSamLineUnique(outputSam->outputFiles[0], &outputSam->samLines[samLineId]);
    }
  }
  else {
    for (unsigned int samLineId = 0; samLineId < outputSam->nSamLines; ++samLineId) {
      ++nReads;
      for (unsigned int fileId = 0; fileId < parameters->nReadsFiles; ++fileId) {
        printSamLineMultiple(outputSam->outputFiles[fileId], &outputSam->samLines[samLineId], fileId, nReads);
      }
    }
  }
  if (pthread_mutex_unlock(outputSam->writeMutex) != 0) {
    fprintf(stderr, "Error, Cannot release write mutex!\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  outputSam->nSamLines = 0;
}

/**
 * Add a new SAM information to the buffer
 */
void addSamLine (outputSam_t *outputSam, unsigned int flag, char *chrName, int64_t pos, bwtint_t nHits, bwtint_t hitId, unsigned int nErrors, char *cigar, char *sequence, char *quality) {
  for (unsigned int fileId = 0; fileId < parameters->nOutputFileNames; ++fileId) {
    copyToSamLine(&outputSam->samLines[outputSam->nSamLines], parameters->nReadsFiles, outputSam->counts, flag, chrName, pos, nHits, hitId, nErrors, cigar, sequence, quality);
  }
  ++outputSam->nSamLines;
  if (outputSam->nSamLines == N_SAM_LINES) {
    writeToSam(outputSam);    
  }
}

/**
 * Add one SAM info (too many hits) to the buffer
 */
void printReadLineManyHits (char *seq, char *qual, bwtint_t nHits, unsigned int nErrors, outputSam_t *outputSam) {
  addSamLine(outputSam, 4, "*", 0, nHits, 0, nErrors, "*", seq, qual);
}

/**
 * Decide whether to use the direct or the reverse-complement before add one SAM info to the buffer
 */
void printReadLine (bool forward, unsigned int flag, char *chrName, int64_t pos, bwtint_t nHits, bwtint_t hitId, unsigned int nErrors, outputSam_t *outputSam) {
  char *cigar, *seq, *qual;
  if (! outputSam->isBackwardSet) {
    computeReverseComplement(outputSam);
  }
  if (forward) {
    cigar = outputSam->forwardCigar;
    seq   = outputSam->forwardSeq;
    qual  = outputSam->forwardQual;
  }
  else {
    cigar = outputSam->backwardCigar;
    seq   = outputSam->backwardSeq;
    qual  = outputSam->backwardQual;
  }
  addSamLine(outputSam, flag, chrName, pos, nHits, hitId, nErrors, cigar, seq, qual);
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

void setCountsSam (outputSam_t *outputSam, count_t *counts) {
  memcpy(outputSam->counts, counts, parameters->nReadsFiles * sizeof(count_t));
}

void openSamFiles (FILE **outputFiles) {
  for (unsigned int fileId = 0; fileId < parameters->nOutputFileNames; ++fileId) {
    outputFiles[fileId] = fopen(parameters->outputSamFileNames[fileId], "w");
    if (outputFiles[fileId] == NULL) {
      fprintf(stderr, "Error, cannot open output file '%s' for writing!\nExiting.\n", parameters->outputSamFileNames[fileId]);
      exit(EXIT_FAILURE);
    }
    fprintf(outputFiles[fileId], "@HD\tVN:1.6\tSO:unsorted\n");
    for (int i = 0; i < bns->n_seqs; ++i) {
      fprintf(outputFiles[fileId], "@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
    }
  }
}

void closeSamFiles (FILE **outputFiles) {
  for (unsigned int fileId = 0; fileId < parameters->nOutputFileNames; ++fileId) {
    fclose(outputFiles[fileId]);
  }
}

#endif
