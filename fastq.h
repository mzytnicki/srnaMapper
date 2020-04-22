#ifndef FASTQ_H
#define FASTQ_H

#include "constants.h"
#include "parameters.h"
#include "tree.h"

 /**
  * Read a FASTQ file, and store it to a tree
  */
void readReadsFile (char *fileName, tree_t *tree, unsigned int fileId) {
  FILE *inFile;
  char *line = NULL;
  char *sequence = NULL;
  char *quality = NULL;
  size_t len = 0;
  ssize_t nRead;
  inFile = fopen(fileName, "r");
  if (inFile == NULL) exit(EXIT_FAILURE);
  while ((nRead = getline(&line, &len, inFile)) != -1) {
    nRead = getline(&sequence, &len, inFile);
    if (nRead == -1) {
      fprintf(stderr, "Input file '%s' is corrupted.\nAborting.\n", fileName);
      exit(EXIT_FAILURE);
    }
    nRead = getline(&line, &len, inFile);
    if (nRead == -1) {
      fprintf(stderr, "Input file '%s' is corrupted.\nAborting.\n", fileName);
      exit(EXIT_FAILURE);
    }
    nRead = getline(&quality, &len, inFile);
    if (nRead == -1) {
      fprintf(stderr, "Input file '%s' is corrupted.\nAborting.\n", fileName);
      exit(EXIT_FAILURE);
    }
    assert(strlen(sequence) == strlen(quality));
    assert(strlen(sequence) == (unsigned long) nRead);
    trimSequence(nRead, sequence);
    trimSequence(nRead, quality);
    ++stats->nReads;
    if (! addSequence(tree, nRead-1, sequence, quality, fileId)) {
      ++stats->nShortReads;
    }
  }
  free(line);
  free(sequence);
  free(quality);
  fclose(inFile);
}

#endif
