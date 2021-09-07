#ifndef READ_NAMES_H
#define READ_NAMES_H

#include <assert.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "constants.h"
#include "helper.h"

/**
 * Read names data stores all the read names into an "efficient" structure.
 * The structure is a set of chars.
 */

/*
typedef struct {
  char       **sets;
  unsigned int setId;
  uint32_t     offset;
} readNames_t;


static const unsigned int SET_SIZE    = 8;
static const unsigned int OFFSET_SIZE = 32 - SET_SIZE;

static const unsigned int N_SETS     = (1 << SET_SIZE)    - 1;
static const unsigned int MAX_OFFSET = (1 << OFFSET_SIZE) - 1;
*/

/**
 * A read name can be retrieved using:
 *  - the id of the set
 *  - offset
 */

/*
typedef struct {
  uint32_t setId:  8;  // SET_SIZE;
  uint32_t offset: 24; // OFFSET_SIZE;
} readNameId_t;

void createReadNames (readNames_t *readNames) {
  readNames->sets    = (char **) mallocOrDie(N_SETS * sizeof(char *));
  readNames->sets[0] = (char *)  mallocOrDie(MAX_OFFSET * sizeof(char));
  readNames->setId   = 0;
  readNames->offset  = 0;
}

void freeReadNames (readNames_t *readNames) {
  for (unsigned int setId = 0; setId < readNames->setId; ++setId) {
    free(readNames->sets[setId]);
  }
  free(readNames->sets);
}

void addReadName (readNames_t *readNames, const char *readName, int readNameLength, readNameId_t *readNameId) {
  // Account for the terminal '\0'
  // printf("Adding name %s of size %i, set %" PRIu32 "/%" PRIu32 ", offset %" PRIu32 "/%" PRIu32 "\n", readName, readNameLength, readNames->setId, N_SETS, readNames->offset, MAX_OFFSET);
  ++readNameLength;
  if (readNames->offset + readNameLength >= MAX_OFFSET) {
    if (readNames->setId == N_SETS) {
      fprintf(stderr, "Cannot store that many read names.  This is a serious problem.  Please contact the developper, matthias.zytnicki@inrae.fr to correct this.\n");
      exit(EXIT_FAILURE);
    }
    ++readNames->setId;
    readNames->sets[readNames->setId] = (char *) mallocOrDie(MAX_OFFSET * sizeof(char));
    readNames->offset = 0;
  }
  memcpy(readNames->sets[readNames->setId] + readNames->offset, readName, readNameLength * (sizeof(char)));
  readNameId->setId  = readNames->setId;
  readNameId->offset = readNames->offset;
  readNames->offset += readNameLength;
}

char* getReadNames (readNames_t *readNames, readNameId_t *readName) {
  return readNames->sets[readName->setId] + readName->offset;
}
*/

/**
 * This concatenates read names, each one being separated by '\0'
 */
typedef struct {
  size_t length;
  char  *readNames;
} readNames_t;

void setReadName (readNames_t *readNames, char *readName, size_t length) {
  // Account for terminal '\0'
  ++length;
  readNames->readNames = (char *) mallocOrDie(length * sizeof(char));
  memcpy(readNames->readNames, readName, length);
  readNames->length = length;
}

void appendReadName (readNames_t *readNames, char *readName, size_t length) {
  // Account for terminal '\0'
  ++length;
  size_t newLength = length + readNames->length;
  readNames->readNames = (char *) reallocOrDie(readNames->readNames, newLength * sizeof(char));
  memcpy(readNames->readNames + readNames->length, readName, length);
  readNames->length = newLength;
}

void addReadName (readNames_t *readNames, char *readName, size_t length) {
  if (readNames->length == 0) {
    setReadName(readNames, readName, length);
  }
  else {
    appendReadName(readNames, readName, length);
  }
}

void freeReadName (readNames_t *readNames) {
  free(readNames->readNames);
}

#endif
