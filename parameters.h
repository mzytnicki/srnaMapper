#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <time.h>
#include "constants.h"

typedef struct {
  bool uniqueOutputFile;
  char **readsFileNames;
  char  *genomeFileName;
  char  *outputReadsFileName;
  char **outputSamFileNames;
  unsigned int nReadsFiles;
  unsigned int nThreads;
  size_t       maxNErrors;
  unsigned int lowComplexityThreshold;
  unsigned int maxNHits;
  unsigned int nOutputFileNames;
} parameters_t;

parameters_t *parameters;

void printUsage () {
  puts("Usage:\n"
       "srnaMapper [parameters]\n\n"
       "Compulsory parameters:\n"
       "  -r string: file name in FASTQ format\n"
       "  -g string: prefix of the genome database (produced by bwa build)\n"
       "  -o string: output file in SAM format\n"
       "Optional parameters:\n"
       "  -e int: maximum number of errors (default: 2)\n"
       "  -t int: number of threads (default: 1)\n"
       "  -n int: discard reads when they map more than n times (default: 5)\n"
       "  -f int: low complexity threshold, more is more lenient (default: 6)\n"
       "  -u: if set, print all the mapped reads in a unique SAM file (with the counts for each sample)\n"
       "  -s int: set the random seed (time otherwise)\n"
       "Notes:\n"
       "  The '-r' option should be repeated once per input file.\n"
       "  Unless the the '-u' option is set, the '-o' option should also be repeated once per input file.\n");
}

void parseCommandLine (int argc, char const **argv) {
  char *endptr;
  srand(time(NULL));
  parameters->nOutputFileNames       = 0;
  parameters->uniqueOutputFile       = false;
  parameters->genomeFileName         = NULL;
  parameters->outputReadsFileName    = NULL;
  parameters->nReadsFiles            = 0;
  parameters->nThreads               = 1;
  parameters->maxNErrors             = 2;
  parameters->lowComplexityThreshold = 6;
  parameters->maxNHits               = 5;
  // Count the number of reads files and output files first
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-r") == 0) {
      ++parameters->nReadsFiles;
    }
    else if (strcmp(argv[i], "-o") == 0) {
      ++parameters->nOutputFileNames;
    }
  }
  parameters->readsFileNames     = (char **) mallocOrDie(parameters->nReadsFiles      * sizeof(char *));
  parameters->outputSamFileNames = (char **) mallocOrDie(parameters->nOutputFileNames * sizeof(char *));
  parameters->nReadsFiles        = 0;
  parameters->nOutputFileNames   = 0;
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-h") == 0) {
      printUsage();
      exit(EXIT_SUCCESS);
    }
    if (strcmp(argv[i], "-v") == 0) {
      printf("%s v%s\n", argv[0], VERSION);
      exit(EXIT_SUCCESS);
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
      parameters->outputSamFileNames[parameters->nOutputFileNames] = strdup(argv[i]);
      ++parameters->nOutputFileNames;
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
    else if (strcmp(argv[i], "-u") == 0) {
      parameters->uniqueOutputFile = true;
    }
    else if (strcmp(argv[i], "-s") == 0) {
      ++i;
      srand(strtol(argv[i], &endptr, 10));
    }
    else {
      printf("Cannot understand parameter '%s'\n", argv[i]);
      printUsage();
      exit(EXIT_FAILURE);
    }
  }
  if (parameters->nReadsFiles == 0) {
    printUsage();
    puts("Reads file is missing.\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  if (parameters->genomeFileName == NULL) {
    printUsage();
    puts("Genome index file is missing.\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  if (parameters->nOutputFileNames == 0) {
    printUsage();
    puts("Output SAM file is missing.\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  if ((parameters->uniqueOutputFile) && (parameters->nOutputFileNames != 1)) {
    printUsage();
    puts("Expect only one output file names.\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  if ((! parameters->uniqueOutputFile) && (parameters->nOutputFileNames != parameters->nReadsFiles)) {
    printUsage();
    printf("The number of input files (%u) and output files (%u) differ.\nExiting.\n", parameters->nReadsFiles, parameters->nOutputFileNames);
    exit(EXIT_FAILURE);
  }
}

void freeParameters(parameters_t *p) {
  for (unsigned int readsFileId = 0; readsFileId < p->nReadsFiles; ++readsFileId) {
    free(p->readsFileNames[readsFileId]);
  }
  if (p->readsFileNames != NULL) {
    free(p->readsFileNames);
  }
  if (p->genomeFileName != NULL) {
    free(p->genomeFileName);
  }
  if (p->outputReadsFileName != NULL) {
    free(p->outputReadsFileName);
  }
  for (unsigned int outputFileId = 0; outputFileId < p->nOutputFileNames; ++outputFileId) {
    free(p->outputSamFileNames[outputFileId]);
  }
  free(p->outputSamFileNames);
}

#endif
