#ifndef STATS_H
#define STATS_H

#include "constants.h"
#include "parameters.h"

typedef struct {
  unsigned long int      nReads;
  unsigned long int      nShortReads;
  //unsigned long int      nDown;
  //unsigned long int      nDownPreprocessed;
  unsigned long int      nBufferCalls;
  unsigned long int      nBufferCallSucesses;
  size_t                 maxHashSize;
  size_t                 maxHashVectorSize;
  size_t                *maxNStates;
  unsigned long int      nTentativeStateInsertions;
  unsigned long int      nSkippedStateInsertions;
  unsigned long long int nBwtPerDepth[MAX_READ_LENGTH];
} stats_t;

stats_t *stats;

void initializeStats () {
  stats->nReads                    = 0;
  stats->nShortReads               = 0;
  stats->nBufferCalls              = 0;
  stats->nBufferCallSucesses       = 0;
  stats->maxHashSize               = 0;
  stats->maxHashVectorSize         = 0;
  stats->maxNStates                = (size_t *) calloc(parameters->maxNErrors+1, sizeof(size_t));
  stats->nTentativeStateInsertions = 0;
  stats->nSkippedStateInsertions   = 0;
  memset(stats->nBwtPerDepth, 0, MAX_READ_LENGTH * sizeof(unsigned long long int));
}

void freeStats () {
  free(stats->maxNStates);
}

void printStats () {
  char *savedLocale;
  savedLocale = setlocale (LC_ALL, NULL);
  setlocale(LC_NUMERIC, "");
  printf("Very small sequences: %'lu/%'lu\n", stats->nShortReads, stats->nReads);
  /*
  printf("# max states:");
  for (size_t nErrors = 0; nErrors <= parameters->maxNErrors; ++nErrors) {
    printf(" %'zu: %'zu ", nErrors, stats->maxNStates[nErrors]);
  }
  printf("\n");
  printf("# buffer call successes %'lu/%'lu (%i%%)\n", stats->nBufferCallSucesses, stats->nBufferCalls, (stats->nBufferCalls == 0)? 0: (int) (round(((double) stats->nBufferCallSucesses) / stats->nBufferCalls * 100)));
  printf("# max state hash: hash %'lu/%'i, vector: %'lu\n", stats->maxHashSize, N_STATES_HASH_SIZE, stats->maxHashVectorSize);
  printf("# state insertion skipped: %'lu/%'lu (%i%%)\n", stats->nSkippedStateInsertions, stats->nTentativeStateInsertions, (stats->nTentativeStateInsertions == 0)? 0: (int) (round(((double) stats->nSkippedStateInsertions) / stats->nTentativeStateInsertions * 100)));
  printf("# BWT calls:\n");
  for (size_t i = 0; i < MAX_READ_LENGTH; ++i) {
    if (stats->nBwtPerDepth[i] > 0) {
      printf("\tDepth %3zu: %'13llu\n", i, stats->nBwtPerDepth[i]);
    }
  }
  */
  setlocale(LC_ALL, savedLocale);
}

#endif
