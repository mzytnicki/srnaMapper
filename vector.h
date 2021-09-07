#ifndef VECTOR_H
#define VECTOR_H

#include "constants.h"
#include "helper.h"

/** 
 * A very simple resizeable generic vector (with holes).
 */
typedef struct {
  void    *data;
  size_t   dataSize;
  uint64_t nElements;
  uint64_t nAllocated;
} vector_t;

void createVector (vector_t *vector, size_t dataSize, uint64_t nElements) {
  vector->dataSize   = dataSize;
  vector->nElements  = nElements;
  vector->nAllocated = nElements;
  vector->data       = callocOrDie(nElements, vector->dataSize);
}

void freeVector (vector_t *vector) {
  free(vector->data);
}

void resizeVector (vector_t *vector, uint64_t index) {
  if (index < vector->nAllocated) {
    return;
  }
  uint64_t previousNAllocated = vector->nAllocated;
  while (index >= vector->nAllocated) {
    vector->nAllocated *= 2;
  }
  vector->data = reallocOrDie(vector->data, vector->nAllocated * vector->dataSize);
  memset(vector->data + previousNAllocated * vector->dataSize, 0, (vector->nAllocated - previousNAllocated) * vector->dataSize);
}

void setVectorElement (vector_t *vector, uint64_t index, void *data) {
  ++vector->nElements;
  //printf("Setting element of size %zu @ %" PRIu64 "/%" PRIu64 "\n", vector->dataSize, index, vector->nAllocated); fflush(stdout);
  resizeVector(vector, index);
  memcpy(vector->data + (index * vector->dataSize), data, vector->dataSize);
}

void copyVectorElement (vector_t *vector, uint64_t index, void *element) {
  if (index < vector->nAllocated) {
    memcpy(vector->data + (index * vector->dataSize), element, vector->dataSize);
  }
}

void *getVectorElement (vector_t *vector, uint64_t index) {
  resizeVector(vector, index);
  return vector->data + (index * vector->dataSize);
}

bool isVectorSet (vector_t *vector, uint64_t index) {
  if (index >= vector->nAllocated) {
    return false;
  }
  return (((char *) vector->data)[index * vector->dataSize] != 0);
}

#endif
