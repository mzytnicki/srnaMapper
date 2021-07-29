#ifndef QUALITY_H
#define QUALITY_H

#include "constants.h"

/******* Quality type *******/
/**
 * A quality the quality of a set of reads with the same sequence.
 * It is:
 *   - the qualities: an array of char*
 *   - the size of the array
 * The quality array is roughly of the same size as the tree array.
 * If the quality is set for a cell, the corresponding index should not be null.
 */

typedef struct {
  char       **qualities;
  //uint64_t    *cellIds;
  uint32_t nQualities;
  uint32_t nAllocated;
} quality_t;

void createQualities (quality_t *qualities) {
  qualities->nAllocated = INIT_N_QUALITIES;
  qualities->qualities  = (char **) callocOrDie(qualities->nAllocated, sizeof(char *));
  qualities->nQualities = 0;
}

void freeQualities (quality_t *qualities) {
  for (uint32_t i = 0; i < qualities->nAllocated; ++i) {
    if (qualities->qualities[i] != NULL) {
      free(qualities->qualities[i]);
    }
  }
  free(qualities->qualities);
}

uint32_t findQualityId (const quality_t *qualities, uint64_t cellId) {
  //printf("\t\t\tFinding quality for %" PRIu64 "\n", cellId); fflush(stdout);
  if ((cellId >= qualities->nAllocated) || (qualities->qualities[cellId] == NULL)) return NO_QUALITY;
  return cellId;
}

char *findQuality (const quality_t *qualities, uint64_t cellId) {
  uint32_t qualityId = findQualityId(qualities, cellId);
  if (qualityId == NO_QUALITY) return NULL;
  //printf("\t\t\tFinding quality for %" PRIu64 ": %u\n", cellId, qualityId); fflush(stdout);
  return qualities->qualities[qualityId]; 
}

void addQuality (quality_t *qualities, uint64_t cellId, size_t l, char *quality) {
  //printf("Adding quality %s at %" PRIu64 ".\n", quality, cellId);
  //TODO: Is the "if" useful?
  size_t previousNAllocated;
  ++qualities->nQualities;
  //printf("\tNow %" PRIu32 " qualities\n", qualities->nQualities);
  if (cellId >= qualities->nAllocated) {
    previousNAllocated = qualities->nAllocated;
    while (cellId >= qualities->nAllocated) {
      qualities->nAllocated *= 2;
    }
    if ((qualities->qualities = (char **) realloc(qualities->qualities, qualities->nAllocated * sizeof(char *))) == NULL) {
      printf("Cannot allocate memory for qualities of size %u.\nExiting.\n", qualities->nAllocated);
      exit(EXIT_FAILURE);
    }
    for (size_t i = previousNAllocated; i < qualities->nAllocated; ++i) {
      qualities->qualities[i] = NULL;
    }
  }
  qualities->qualities[cellId] = strndup(quality, l);
}

void mergeQualities (char *quality1, const char *quality2, size_t l) {
  assert(strlen(quality1) == strlen(quality2));
  assert(strlen(quality1) == l);
  for (size_t i = 0; i < l; ++i) {
    quality1[i] = MAX(quality1[i], quality2[i]);
  }
}

void replaceQuality (quality_t *qualities, uint32_t qualityId, size_t l, char *quality) {
  //printf("Quality: %zu vs %zu\n", strlen(qualities->qualities[qualityId]), strlen(quality));
  assert(strlen(quality) == l);
  mergeQualities(qualities->qualities[qualityId], quality, l);
}

void _setQuality (quality_t *qualities, uint64_t cellId, size_t l, char *quality) {
  uint32_t qualityId = findQualityId(qualities, cellId);
  if (qualityId == NO_QUALITY) {
    addQuality(qualities, cellId, l, quality);
  }
  else {
    replaceQuality(qualities, qualityId, l, quality);
  }
}


#endif
