#ifndef CELL_INFO_H
#define CELL_INFO_H

#include <assert.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "constants.h"

#define NO_INFO ((uint32_t) -1)

/**
 * CellInfo stores informations on compressed tree cells:
 *  - the id of the cell
 *  - the quality
 *  - the counts
 */
/*
typedef struct {
  uint32_t      cellId;
  char         *quality;
  count_t      *counts;
  char        **readNames;
} cellInfo_t;

void createEmptyCellInfo (cellInfo_t *cellInfo) {
  cellInfo->cellId = NO_INFO;
}

bool isEmptyCellInfo (const cellInfo_t *cellInfo) {
  return (cellInfo->cellId == NO_INFO);
}

void createCellInfo (cellInfo_t *cellInfo, size_t readSize, size_t nSamples) {
  ++readSize;
  cellInfo->quality   = (char *)    mallocOrDie(readSize * sizeof(char));
  cellInfo->counts    = (count_t *) mallocOrDie(nSamples * sizeof(count_t));
  cellInfo->readNames = (char **)   mallocOrDie(nSamples * sizeof(char *));
}

void freeCellInfo (cellInfo_t *cellInfo) {
  free(cellInfo->quality);
  free(cellInfo->counts);
  for (unsigned int readFileId = 0; readFileId < parameters->nReadsFiles; ++readFileId) {
    free(cellInfo->readNames[readFileId]);
  }
  free(cellInfo->readNames);
}

void setCellInfo (cellInfo_t *cellInfo, uint32_t cellId, char *quality, size_t readSize, count_t *counts, unsigned int nSamples) {
  cellInfo->cellId = cellId;
  ++readSize;
  memcpy(cellInfo->quality, quality, readSize * sizeof(char));
  memcpy(cellInfo->counts,  counts,  nSamples * sizeof(count_t));
}

void setCellInfoReadNames (cellInfo_t *cellInfo, char *readNames, size_t length, unsigned int sampleId) {
  strndup(cellInfo->readNames[sampleId], readNames, length);
}

void printCellInfo (const cellInfo_t *cellInfo) {
  printf("cell id: %" PRIu32 ", qual: %s, counts (%p): ", cellInfo->cellId, cellInfo->quality, cellInfo->counts);
  for (unsigned int countId = 0; countId < parameters->nReadsFiles; ++countId) {
    printf("%" PRIu32 " ", cellInfo->counts[countId]);
  }
  printf("\n");
}
*/

/**
 * CellInfos store all the cell infos in a vector.
 * The infos a supposed to be read sequentially, from the first cell to the last (possibly skipping some)...
 */
typedef struct {
  size_t        nAllocatedCellInfos;
  size_t        nCellInfos;
  uint32_t     *cellIds;
  char        **qualities;
  count_t      *counts;
  char        **readNames;
} cellInfos_t;

typedef size_t cellVisitor_t;

void clearCellVisitor (cellVisitor_t *cellVisitor) {
  *cellVisitor = 0;
}

void createCellInfos (cellInfos_t *cellInfos, size_t nInfos /* , size_t readSize, size_t nSamples */) {
  cellInfos->nAllocatedCellInfos = nInfos;
  cellInfos->nCellInfos          = 0;
  cellInfos->cellIds             = (uint32_t *) mallocOrDie((nInfos+1) * sizeof(uint32_t));  
  cellInfos->counts              = (count_t *) mallocOrDie(nInfos * (parameters->nReadsFiles) * sizeof(count_t));  
  cellInfos->qualities           = (char **) mallocOrDie(nInfos * sizeof(char *));  
  cellInfos->readNames           = (char **) mallocOrDie(nInfos * sizeof(char *));  
  cellInfos->cellIds[nInfos]     = NO_INFO;
}

void freeCellInfos (cellInfos_t *cellInfos) {
  for (size_t cellInfosId = 0; cellInfosId < cellInfos->nCellInfos; ++cellInfosId) {
    free(cellInfos->qualities[cellInfosId]);
  }
  free(cellInfos->cellIds);
  free(cellInfos->qualities);
  free(cellInfos->counts);
  for (size_t cellInfosId = 0; cellInfosId < cellInfos->nCellInfos; ++cellInfosId) {
    free(cellInfos->readNames[cellInfosId]);
  }
  free(cellInfos->readNames);
}

/**
 * Allocate memory for a new cellInfo, and copy content.
 * Returns the buffer of the read names.
 */
char *addCellInfo (cellInfos_t *cellInfos, uint32_t cellId, char *quality, size_t readSize, count_t *counts, size_t readNameLength) {
  assert(cellInfos->nCellInfos < cellInfos->nAllocatedCellInfos);
  //printf("Setting cell info #%" PRIu32 "\n", cellId);
  cellInfos->cellIds[cellInfos->nCellInfos] = cellId;
  memcpy(cellInfos->counts + (parameters->nReadsFiles * cellInfos->nCellInfos), counts, parameters->nReadsFiles * sizeof(count_t));
  cellInfos->qualities[cellInfos->nCellInfos] = (char *) malloc(readSize + 1);
  if (cellInfos->qualities[cellInfos->nCellInfos] == NULL) {
    fprintf(stderr, "Error, cannot add quality to compressed tree.\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  memcpy(cellInfos->qualities[cellInfos->nCellInfos], quality, readSize+1);
  cellInfos->readNames[cellInfos->nCellInfos] = (char *) mallocOrDie(readNameLength);
  ++cellInfos->nCellInfos;
  return cellInfos->readNames[cellInfos->nCellInfos-1];
}

//cellInfo_t *getCellInfo (uint32_t cellId, cellVisitor_t *cellVisitor) {
bool getCellInfo (const cellInfos_t *cellInfos, uint32_t cellId, cellVisitor_t *cellVisitor) {
  //printf("Looking for cell info %" PRIu32 ", with visitor %zu and %zu cells\n", cellId, *cellVisitor, cellInfos->nCellInfos); fflush(stdout);
  //printf("Current state is %zu -> %" PRIu32 "\n", *cellVisitor, cellInfos->cellIds[*cellVisitor]); fflush(stdout);
  /*
  for (; (! isEmptyCellInfo(cellVisitor->currentCellInfo)) && (cellVisitor->currentCellInfo->cellId < cellId); ++cellVisitor->currentCellInfo) {
    //printf("Got %" PRIu32 "\n", cellVisitor->currentCellInfo->cellId); fflush(stdout);
  }
  return (cellVisitor->currentCellInfo->cellId == cellId)? cellVisitor->currentCellInfo: NULL;
  */
  for (; cellInfos->cellIds[*cellVisitor] < cellId; ++(*cellVisitor)) {
    if (cellInfos->cellIds[*cellVisitor] == NO_INFO) {
      //printf("Over\n"); fflush(stdout);
      return false;
    }
    //printf("Got %" PRIu32 "\n", cellInfos->cellIds[*cellVisitor]); fflush(stdout);
  }
  return (cellInfos->cellIds[*cellVisitor] == cellId);
}

count_t *getCounts (const cellInfos_t *cellInfos, cellVisitor_t cellVisitor) {
  return cellInfos->counts + (cellVisitor * parameters->nReadsFiles);
}

char *getQuality (const cellInfos_t *cellInfos, cellVisitor_t cellVisitor) {
  return cellInfos->qualities[cellVisitor];
}

char *getCellInfoReadNames (const cellInfos_t *cellInfos, cellVisitor_t cellVisitor) {
  return cellInfos->readNames[cellVisitor];
}

#endif
