#ifndef FASTQ_H
#define FASTQ_H

#include "constants.h"
#include "parameters.h"
#include "tree.h"

#define FASTQ_LINES 4

/**
 * Open/check the input FASTQ files)
 */
void openFastqFiles (FILE **inputFastqFiles) {
  for (size_t fileId = 0; fileId < parameters->nReadsFiles; ++fileId) {
    inputFastqFiles[fileId] = fopen(parameters->readsFileNames[fileId], "r");
    if (inputFastqFiles[fileId] == NULL) {
      fprintf(stderr, "Error, cannot open FASTQ file '%s'!\nExiting.\n", parameters->readsFileNames[fileId]);
      exit(EXIT_FAILURE);
    }
  }
}

void closeFastqFiles (FILE **inputFastqFiles) {
  for (size_t fileId = 0; fileId < parameters->nReadsFiles; ++fileId) {
    fclose(inputFastqFiles[fileId]);
  }
}

/**
 * Jump to the next FASTQ chunk.
 * Return false if we reached the EOF.
 * Should read the first 5 lines, because the first one could be truncated.
 */
bool setToStart (FILE *inFile, off64_t posStart) {
  char *line = NULL;
  size_t len = 0;
  char    firstChars[2*FASTQ_LINES];
  off64_t offsets[2*FASTQ_LINES-1];
  // Bad trick in case the 5th line is over the EOF.
  firstChars[2*FASTQ_LINES-1] = 0;
  if (fseeko64(inFile, posStart, SEEK_SET) != 0) {
    fprintf(stderr, "Error, problem reading FASTQ file.\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  for (unsigned int lineId = 0; lineId < 2*FASTQ_LINES-1; ++lineId) {
    offsets[lineId] = ftello64(inFile);
    if (getline(&line, &len, inFile) == -1) {
      // If the 5th line does not exist, it might still be ok...
      if (lineId == FASTQ_LINES) {
        break;
      }
      if (line != NULL) free(line);
      return false;
    }
    firstChars[lineId] = line[0];
    if (lineId >= 2) {
      if ((firstChars[lineId] == '+') && (firstChars[lineId-2] == '@')) {
        if (fseeko64(inFile, offsets[lineId-2], SEEK_SET) != 0) {
          fprintf(stderr, "Error, problem reading FASTQ file.\nExiting.\n");
          exit(EXIT_FAILURE);
        }
        if (line != NULL) free(line);
        return true;
      }
    }
  }
  fprintf(stderr, "Error, cannot find the beginning of the chunk of the FASTQ file.\nExiting.\n");
  if (line != NULL) free(line);
  exit(EXIT_FAILURE);
  return false;
}

bool isPastEnd (off64_t posCurrent, off64_t posEnd) {
  return (posCurrent > posEnd);
}


/**
 * Read a FASTQ file, and store it to a tree
 * Return false if we reach EOF
 */
bool readReadsFile (FILE *inFile, char *fileName, tree_t *tree, unsigned int fileId, off64_t posStart, off64_t posEnd) {
  char *line = NULL;
  char *sequence = NULL;
  char *quality = NULL;
  size_t len = 0;
  bool over = false;
  ssize_t nRead;
  if (! setToStart(inFile, posStart)) {
    return false;
  }
  while (((nRead = getline(&line, &len, inFile)) != -1) && (! over)) {
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
    addSequence(tree, nRead-1, sequence, quality, fileId);
    over = isPastEnd(ftello64(inFile), posEnd);
  }
  free(line);
  free(sequence);
  free(quality);
  return over;
}

#endif
