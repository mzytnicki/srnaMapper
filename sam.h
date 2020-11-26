#ifndef SAM_H
#define SAM_H

#include "constants.h"
#include "helper.h"
#include "parameters.h"

/**
 * Stores a hit, and the number of corresponding identical sequences in the input files
 */
typedef struct {
  char         threadId; // one char ('A', 'B', etc.) per thread, written in the (unique) read name
  count_t     *readIds;  // one read id for each output file, the id is given to the first read (in case of multiple reads per sequence)
  count_t     *counts;
  unsigned int flag;
  int          rid;
  uint64_t     pos;
  char         cigar[MAX_CIGAR_SIZE];
  char        *sequence;
  char        *quality;
  unsigned int nHits;
  unsigned int hitId;
  unsigned int nErrors;
} samLine_t;

void createSamLine (samLine_t *samLine, char threaId, count_t *readIds, count_t *counts, char *sequence, char *quality) {
  samLine->threadId = threaId;
  samLine->readIds  = readIds;
  samLine->counts   = counts;
  samLine->sequence = sequence;
  samLine->quality  = quality;
}

/*
void freeSamLine (samLine_t *samLine) {
  free(samLine->counts);
  free(samLine->sequence);
  free(samLine->quality);
}
*/

int samLineComparator (const void *samLine1, const void *samLine2) {
  int difference = ((samLine_t *) samLine1)->rid - ((samLine_t *) samLine2)->rid;
  if (difference != 0) return difference;
  return ((samLine_t *) samLine1)->pos - ((samLine_t *) samLine2)->pos;
}

unsigned int computeMapq (unsigned int nHits, unsigned int nErrors) {
  if ((nHits > 1) || (nErrors >= 40)) {
    return 0;
  }
  else {
    return 40 - nErrors;
  }
}

void swapSamLines (samLine_t *samLine1, samLine_t *samLine2) {
  samLine_t tmp = *samLine1;
  *samLine1 = *samLine2;
  *samLine2 = tmp;
}

void printSamLine (const samLine_t *samLine) {
  printf("%s (%p, %p) pos %i:%" PRId64 " %u/%u hits, counts:", samLine->sequence, samLine, samLine->sequence, samLine->rid, samLine->pos, samLine->hitId, samLine->nHits);
  for (unsigned int i = 0; i < parameters->nReadsFiles; ++i) {
    printf(" %u", samLine->counts[i]);
  }
  printf("\n");
}

/**
 * Copy one SAM information to the buffer
 */
void copyToSamLine (samLine_t *samLine, size_t nCounts, count_t *readIds, count_t *counts, unsigned int flag, int rid, int64_t pos, bwtint_t nHits, bwtint_t hitId, unsigned int nErrors, char *cigar, char *sequence, char *quality) {
  memcpy(samLine->readIds,  readIds,  nCounts * sizeof(count_t));
  memcpy(samLine->counts,   counts,   nCounts * sizeof(count_t));
  strcpy(samLine->cigar,    cigar);
  strcpy(samLine->sequence, sequence);
  strcpy(samLine->quality,  quality);
  samLine->rid     = rid;
  samLine->flag    = flag;
  samLine->pos     = pos;
  samLine->nHits   = nHits;
  samLine->hitId   = hitId;
  samLine->nErrors = nErrors;
  //printf("Adding sam line: ");
  //printSamLine(samLine);
} 

void printSamTailLine (FILE *file, samLine_t *samLine) {
  fprintf(file, "\t%u\t%s\t%" PRId64 "\t%d\t%s\t*\t0\t0\t%s\t%s\tNH:i:%u\tHI:i:%u\tIH:i:%u\tNM:i:%u\n", samLine->flag, (samLine->rid == -1)? "*": bns->anns[samLine->rid].name, samLine->pos, computeMapq(samLine->nHits, samLine->nErrors), samLine->cigar, samLine->sequence, samLine->quality, samLine->nHits, samLine->hitId, samLine->nHits, samLine->nErrors);
}

void printSamLineUnique (FILE *file, samLine_t *samLine) {
  fprintf(file, "read_%c_%u_x", samLine->threadId, samLine->readIds[0]);
  for (size_t readsFileId = 0; readsFileId < parameters->nReadsFiles; ++readsFileId) {
    fprintf(file, "_%" PRIu32, samLine->counts[readsFileId]);
  }
  printSamTailLine(file, samLine);
}

void printSamLineMultiple (FILE *file, samLine_t *samLine, unsigned int fileId) {
  for (unsigned int countId = 0; countId < samLine->counts[fileId]; ++countId) {
    fprintf(file, "read_%c_%u", samLine->threadId, samLine->readIds[fileId] + countId);
    printSamTailLine(file, samLine);
  }
}

typedef struct {
  count_t          *readIds;
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
  size_t           maxSamLines;
  pthread_mutex_t *writeMutex;
  samLine_t       *samLines;
  count_t         *samLinesReadIds;
  count_t         *samLinesCounts;
  char            *samLinesSequences;
  char            *samLinesQualities;
  size_t           nSamLines;
} outputSam_t;

void createOutputSam (outputSam_t *outputSam, unsigned int readSize, pthread_mutex_t *mutex, char threadId) {
  count_t *samLinesReadIds;
  count_t *samLinesCounts;
  char    *samLinesSequences;
  char    *samLinesQualities;
  unsigned int readIdOffset    = (parameters->uniqueOutputFile)? 1: parameters->nReadsFiles;
  outputSam->maxSamLines       = MAX(N_SAM_LINES, 2 * parameters->maxNHits);
  outputSam->sequenceSize      = readSize;
  outputSam->isBackwardSet     = false;
  outputSam->readIds           = (count_t *)       calloc(parameters->nReadsFiles,  sizeof(count_t));
  outputSam->counts            = (count_t *)       malloc(parameters->nReadsFiles * sizeof(count_t));
  outputSam->backtraceCigar    = (char *)          malloc(MAX_CIGAR_SIZE          * sizeof(char));
  outputSam->backtraceLengths  = (unsigned char *) malloc(MAX_CIGAR_SIZE          * sizeof(unsigned char));
  outputSam->forwardCigar      = (char *)          malloc(MAX_CIGAR_SIZE          * sizeof(char));
  outputSam->backwardCigar     = (char *)          malloc(MAX_CIGAR_SIZE          * sizeof(char));
  outputSam->forwardSeq        = (char *)          malloc((readSize+1)            * sizeof(char));
  outputSam->backwardSeq       = (char *)          malloc((readSize+1)            * sizeof(char));
  outputSam->forwardQual       = (char *)          malloc((readSize+1)            * sizeof(char));
  outputSam->backwardQual      = (char *)          malloc((readSize+1)            * sizeof(char));
  outputSam->samLinesReadIds   = (count_t *)       malloc(outputSam->maxSamLines * readIdOffset * sizeof(count_t));
  outputSam->samLinesCounts    = (count_t *)       malloc(outputSam->maxSamLines * parameters->nReadsFiles * sizeof(count_t));
  outputSam->samLinesSequences = (char *)          malloc(outputSam->maxSamLines * (readSize+1) * sizeof(char));
  outputSam->samLinesQualities = (char *)          malloc(outputSam->maxSamLines * (readSize+1) * sizeof(char));
  outputSam->writeMutex        = mutex;
  outputSam->nSamLines         = 0;
  outputSam->samLines          = (samLine_t *) malloc(outputSam->maxSamLines * sizeof(samLine_t));
  samLinesReadIds   = outputSam->samLinesReadIds;
  samLinesCounts    = outputSam->samLinesCounts;
  samLinesSequences = outputSam->samLinesSequences;
  samLinesQualities = outputSam->samLinesQualities;
  for (size_t samLineId = 0; samLineId < outputSam->maxSamLines; ++samLineId) {
    createSamLine(&outputSam->samLines[samLineId], threadId, samLinesReadIds, samLinesCounts, samLinesSequences, samLinesQualities);
    samLinesReadIds   += readIdOffset;
    samLinesCounts    += parameters->nReadsFiles;
    samLinesSequences += readSize+1;
    samLinesQualities += readSize+1;
  }
}

void freeOutputSam (outputSam_t *outputSam) {
  free(outputSam->readIds);
  free(outputSam->counts);
  free(outputSam->backtraceCigar);
  free(outputSam->backtraceLengths);
  free(outputSam->forwardCigar);
  free(outputSam->backwardCigar);
  free(outputSam->forwardSeq);
  free(outputSam->backwardSeq);
  free(outputSam->forwardQual);
  free(outputSam->backwardQual);
  free(outputSam->samLinesReadIds);
  free(outputSam->samLinesCounts);
  free(outputSam->samLinesSequences);
  free(outputSam->samLinesQualities);
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
  sprintf(outputSam->forwardCigar, "%u=", readLength);
  sprintf(outputSam->backwardCigar, "%u=", readLength);
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
void writeToSam (outputSam_t *outputSam, bool compulsory) {
  if ((! compulsory) && (outputSam->nSamLines < outputSam->maxSamLines / 2)) {
    return;
  }
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
      for (unsigned int fileId = 0; fileId < parameters->nReadsFiles; ++fileId) {
        //printf("printing output: ");
        //printSamLine(&outputSam->samLines[samLineId]);
        printSamLineMultiple(outputSam->outputFiles[fileId], &outputSam->samLines[samLineId], fileId);
      }
    }
  }
  if (pthread_mutex_unlock(outputSam->writeMutex) != 0) {
    fprintf(stderr, "Error, Cannot release write mutex!\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  outputSam->nSamLines = 0;
}

void removeDuplicatesOutputLines(outputSam_t *outputSam, size_t nLines) {
  assert(nLines <= outputSam->nSamLines);
  size_t       samLineOffset  = outputSam->nSamLines - nLines;
  size_t       nNewLines      = 1;
  samLine_t   *samLineFirst   = outputSam->samLines + samLineOffset;
  samLine_t   *samLineFree    = samLineFirst + 1;
  samLine_t   *samLineCurrent = samLineFree;
  int          previousRid, currentRid;
  uint64_t     previousPos, currentPos;
  unsigned int previousStrand, currentStrand;
  if (nLines <= 1) {
    return;
  }
  /*
  printf("Reducing\n");
  for (unsigned int i = 0; i < nLines; ++i) {
    printSamLine(&outputSam->samLines[samLineOffset+i]);
  }
  */
  qsort(samLineFirst, nLines, sizeof(samLine_t), samLineComparator);
  previousRid    = samLineFirst->rid;
  previousPos    = samLineFirst->pos;
  previousStrand = samLineFirst->flag & CIGAR_REVERSE;
  /*
  printf("Sorting\n");
  for (unsigned int i = 0; i < nLines; ++i) {
    printSamLine(&outputSam->samLines[samLineOffset+i]);
  }
  */
  for (unsigned int samLineId = 1; samLineId < nLines; ++samLineId, ++samLineCurrent) {
    /*
    printf("Step: previousRid: %i, previousPos: %" PRId64 ", current: %p, free: %p\n", previousRid, previousPos, samLineCurrent, samLineFree);
    for (unsigned int i = 0; i < nNewLines; ++i) {
      printSamLine(&outputSam->samLines[samLineOffset+i]);
    }
    */
    currentRid    = samLineCurrent->rid;
    currentPos    = samLineCurrent->pos;
    currentStrand = samLineCurrent->flag & CIGAR_REVERSE;
    if ((previousRid != currentRid) || (currentPos - previousPos > samLineCurrent->nErrors) || (previousStrand != currentStrand)) {
      if (samLineFree != samLineCurrent) {
        //memcpy(samLineFree, samLineCurrent, sizeof(samLine_t));
        swapSamLines(samLineFree, samLineCurrent);
      }
      ++samLineFree;
      ++nNewLines;
      previousRid    = currentRid;
      previousPos    = currentPos;
      previousStrand = currentStrand;
    }
  }
  if (nNewLines != nLines) {
    for (unsigned int hitId = 0; hitId < nNewLines; ++hitId) {
      outputSam->samLines[samLineOffset+hitId].hitId = hitId + 1;
      outputSam->samLines[samLineOffset+hitId].nHits = nNewLines;
    }
  }
  /*
  printf("  to:\n");
  for (unsigned int i = 0; i < nNewLines; ++i) {
    printSamLine(&outputSam->samLines[samLineOffset+i]);
  }
  */
  outputSam->nSamLines = samLineOffset + nNewLines;
  writeToSam(outputSam, false);
}

/**
 * Add a new SAM information to the buffer
 */
void addSamLine (outputSam_t *outputSam, unsigned int flag, int rid, int64_t pos, bwtint_t nHits, bwtint_t hitId, unsigned int nErrors, char *cigar, char *sequence, char *quality) {
  assert(outputSam->nSamLines < outputSam->maxSamLines);
  for (unsigned int fileId = 0; fileId < parameters->nOutputFileNames; ++fileId) {
    copyToSamLine(&outputSam->samLines[outputSam->nSamLines], parameters->nReadsFiles, outputSam->readIds, outputSam->counts, flag, rid, pos, nHits, hitId, nErrors, cigar, sequence, quality);
  }
  ++outputSam->nSamLines;
}

/**
 * Add one SAM info (too many hits) to the buffer
 */
void printReadLineManyHits (char *seq, char *qual, bwtint_t nHits, unsigned int nErrors, outputSam_t *outputSam) {
  addSamLine(outputSam, 4, -1, 0, nHits, 0, nErrors, "*", seq, qual);
  writeToSam(outputSam, false);
  for (unsigned int fileId = 0; fileId < parameters->nOutputFileNames; ++fileId) {
    outputSam->readIds[fileId] += outputSam->counts[fileId];
  }
}

/**
 * Decide whether to use the direct or the reverse-complement before add one SAM info to the buffer
 */
void printReadLine (bool forward, unsigned int flag, int rid, int64_t pos, bwtint_t nHits, bwtint_t hitId, unsigned int nErrors, outputSam_t *outputSam) {
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
  addSamLine(outputSam, flag, rid, pos, nHits, hitId, nErrors, cigar, seq, qual);
}

/**
 * Print a read, which maps only once, and with no error
 */
void printReadUniqueNoError (int strand, int rid, int64_t pos, size_t readLength, outputSam_t *outputSam) {
  unsigned int flag = 0;
  setCigarNoError(outputSam, readLength);
  if (strand == 0) {
    flag = CIGAR_REVERSE;
  }
  printReadLine(strand != 0, flag, rid, pos, 1, 1, 0, outputSam);
  //fprintf(outputSamFile, "%s\t%u\t%s\t%" PRId64 "\t40\t%zuM\t*\t0\t0\t%s\t%s\tNH:i:1\tHI:i:1\tIH:i:1\tNM:i:0\n", qname, flag, chrName, pos, readLength, seq, qual);
  writeToSam(outputSam, false);
}

/**
 * Collect information before printing read
 */
void preparePrintRead (uint64_t pos, int rid, int strand, bwtint_t nHits, bwtint_t hitId, unsigned int nErrors, outputSam_t *outputSam) {
  unsigned int flag = (hitId == 0)? 0: CIGAR_SECONDARY_HIT;
  if (strand == 0) {
    flag |= CIGAR_REVERSE;
  }
  printReadLine(strand != 0, flag, rid, pos, nHits, hitId, nErrors, outputSam);
}

void setCountsSam (outputSam_t *outputSam, count_t *counts) {
  memcpy(outputSam->counts, counts, parameters->nReadsFiles * sizeof(count_t));
}

void updateReadIds (outputSam_t *outputSam) {
  if (parameters->uniqueOutputFile) {
    ++outputSam->readIds[0];
    return;
  }
  for (unsigned int fileId = 0; fileId < parameters->nOutputFileNames; ++fileId) {
    outputSam->readIds[fileId] += outputSam->counts[fileId];
  }
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
