#ifndef STATS_H
#define STATS_H

#include "constants.h"

typedef struct {
  unsigned long int nReads;
  unsigned long int nShortReads;
  unsigned long int nDown;
  unsigned long int nDownPreprocessed;
  unsigned long int nBufferCalls;
  unsigned long int nBufferCallSucesses;
  size_t            maxNStates;
  unsigned long long int nBwtPerDepth[MAX_READ_LENGTH];
} stats_t;

stats_t *stats;

void initializeStats () {
  stats->nReads              = 0;
  stats->nShortReads         = 0;
  stats->nDown               = 0;
  stats->nDownPreprocessed   = 0;
  stats->nBufferCalls        = 0;
  stats->nBufferCallSucesses = 0;
  stats->maxNStates          = 0;
  memset(stats->nBwtPerDepth, 0, MAX_READ_LENGTH * sizeof(unsigned long long int));
}

void printStats () {
  char *savedLocale;
  savedLocale = setlocale (LC_ALL, NULL);
  setlocale(LC_NUMERIC, "");
  printf("Very small sequences: %'lu/%'lu\n", stats->nShortReads, stats->nReads);
  printf("# max states %'zu/%'i\n", stats->maxNStates, N_STATES);
  printf("# buffer call successes %'lu/%'lu (%i%%)\n", stats->nBufferCallSucesses, stats->nBufferCalls, (stats->nBufferCalls == 0)? 0: (int) (round(((double) stats->nBufferCallSucesses) / stats->nBufferCalls * 100)));
  printf("# BWT calls:\n");
  for (size_t i = 0; i < MAX_READ_LENGTH; ++i) {
    if (stats->nBwtPerDepth[i] > 0) {
      printf("\tDepth %3zu: %llu\n", i, stats->nBwtPerDepth[i]);
    }
  }
  setlocale(LC_ALL, savedLocale);
}

#endif
