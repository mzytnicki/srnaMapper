#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "constants.h"

typedef struct {
    char *readsFileNames[255];
    char *genomeFileName;
    char *outputReadsFileName;
    char *outputSamFileName;
    unsigned int nReadsFiles;
    unsigned int nThreads;
    size_t       maxNErrors;
    unsigned int lowComplexityThreshold;
    unsigned int maxNHits;
} parameters_t;

parameters_t *parameters;

void printUsage () {
  puts("srnaCollapser [-h] -r reads -g genome -o filename [-t #threads] [-c filename] [-f filter] [-e #errors] [-n #max_hits]");
}

int parseCommandLine (int argc, char const **argv) {
  char *endptr;
  parameters->genomeFileName         = NULL;
  parameters->outputReadsFileName    = NULL;
  parameters->outputSamFileName      = NULL;
  parameters->nReadsFiles            = 0;
  parameters->nThreads               = 1;
  parameters->maxNErrors             = 2;
  parameters->lowComplexityThreshold = 6;
  parameters->maxNHits               = 5;
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-h") == 0) {
      printUsage();
      return EXIT_FAILURE;
    }
    if (strcmp(argv[i], "-r") == 0) {
      ++i;
      parameters->readsFileNames[parameters->nReadsFiles] = strdup(argv[i]);
      ++parameters->nReadsFiles;
    }
    else if (strcmp(argv[i], "-g") == 0) {
      ++i;
      parameters->genomeFileName = strdup(argv[i]);
    }
    else if (strcmp(argv[i], "-o") == 0) {
      ++i;
      parameters->outputSamFileName = strdup(argv[i]);
    }
    else if (strcmp(argv[i], "-c") == 0) {
      ++i;
      parameters->outputReadsFileName = strdup(argv[i]);
    }
    else if (strcmp(argv[i], "-f") == 0) {
      ++i;
      parameters->lowComplexityThreshold = strtol(argv[i], &endptr, 10);
    }
    else if (strcmp(argv[i], "-e") == 0) {
      ++i;
      parameters->maxNErrors = strtol(argv[i], &endptr, 10);
    }
    else if (strcmp(argv[i], "-t") == 0) {
      ++i;
      parameters->nThreads = strtol(argv[i], &endptr, 10);
    }
    else if (strcmp(argv[i], "-n") == 0) {
      ++i;
      parameters->maxNHits = strtol(argv[i], &endptr, 10);
    }
    else {
      printf("Cannot understand parameter '%s'\n", argv[i]);
      printUsage();
      return EXIT_FAILURE;
    }
  }
  if (parameters->nReadsFiles == 0) {
    printUsage();
    puts("Reads file is missing.\nExiting.");
    return EXIT_FAILURE;
  }
  if (parameters->genomeFileName == NULL) {
    printUsage();
    puts("Genome index file is missing.\nExiting.");
    return EXIT_FAILURE;
  }
  if (parameters->outputSamFileName == NULL) {
    printUsage();
    puts("Output SAM file is missing.\nExiting.");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

void freeParameters(parameters_t *p) {
  for (unsigned int readsFileId = 0; readsFileId < p->nReadsFiles; ++readsFileId) {
    free(p->readsFileNames[readsFileId]);
  }
  if (p->genomeFileName != NULL) {
    free(p->genomeFileName);
  }
  if (p->outputReadsFileName != NULL) {
    free(p->outputReadsFileName);
  }
  if (p->outputSamFileName != NULL) {
    free(p->outputSamFileName);
  }
}

#endif
