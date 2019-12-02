#ifndef TREE2_H
#define TREE2_H

#include <assert.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "quality.h"

#define NO_INFO ((uint32_t) -1)

typedef struct {
  char *read;
  char *quality;
  unsigned int *counts;
} cellInfo_t;

void createCellInfo (cellInfo_t *cellInfo, size_t readSize, size_t nSamples) {
  ++readSize;
  cellInfo->read    = (char *) malloc(readSize * sizeof(char));
  cellInfo->quality = (char *) malloc(readSize * sizeof(char));
  cellInfo->counts  = (unsigned int *) malloc(nSamples * sizeof(unsigned int));
}

void freeCellInfo (cellInfo_t *cellInfo) {
  free(cellInfo->read);
  free(cellInfo->quality);
  free(cellInfo->counts);
}

void setCellInfo (cellInfo_t *cellInfo, char *read, char *quality, size_t readSize, unsigned int *counts, unsigned int nSamples) {
  ++readSize;
  memcpy(&cellInfo->read,    read,    readSize * sizeof(char));
  memcpy(&cellInfo->quality, quality, readSize * sizeof(char));
  memcpy(&cellInfo->counts,  counts,  nSamples * sizeof(unsigned int));
}


typedef struct {
  size_t nAllocatedCellInfos;
  size_t nCellInfos;
  size_t cellInfoId;
  cellInfo_t *cellInfos;
} cellInfos_t;

void createCellInfos (cellInfos_t *cellInfos, size_t nInfos, size_t readSize, size_t nSamples) {
  cellInfos->nAllocatedCellInfos = nInfos;
  cellInfos->nCellInfos = 0;
  cellInfos->cellInfoId = 0;
  cellInfos->cellInfos = (cellInfo_t *) malloc(nInfos * sizeof(cellInfo_t));  
  for (unsigned int i = 0; i < nInfos; ++i) {
    createCellInfo(&cellInfos->cellInfos[i], readSize, nSamples);
  }
}

void freeCellInfos (cellInfos_t *cellInfos) {
  for (unsigned int i = 0; i < cellInfos->nCellInfos; ++i) {
    free(&cellInfos->cellInfos[i]);
  }
  free(cellInfos->cellInfos);
}

void addCellInfo (cellInfos_t *cellInfos, char *read, char *quality, size_t readSize, unsigned int *counts, unsigned int nSamples) {
  assert(cellInfos->nCellInfos < cellInfos->nAllocatedCellInfos);
  setCellInfo(&cellInfos->cellInfos[cellInfos->nCellInfos], read, quality, readSize, counts, nSamples);
  ++cellInfos->nCellInfos;
}

uint32_t getCellInfoId (cellInfos_t *cellInfos, char *read) {
  while(true) {
    if (cellInfos->cellInfoId == cellInfos->nCellInfos) {
      return NO_INFO;
    }
    if (strcmp(read, &cellInfos->cellInfos[cellInfos->cellInfoId].read) != 0) {
      return cellInfos->cellInfoId;
    }
    ++cellInfos->cellInfoId;
  }
}


typedef struct {
  uint32_t nucleotide: 2;
  uint32_t sequence:   26;
  uint32_t length:     4;
  uint32_t cellId;    
} edge2_t;

typedef struct {
  uint32_t      firstEdge;
  unsigned char nEdges;
} cell2_t;


/******* Compressed tree type *******/
/**
 * A tree is stores all the prefixes of the reads.
 * It is:
 *   - a depth
 *   - the number of cells
 *   - the size of the array
 *   - an instance of the quality structure
 * The first N_TREE_BASE cells correspond to the all the combinations of TREE_BASE_SIZE-mers.
 * In other words, the 1-mers, 2-mers, ..., TREE_BASE_SIZE-1-mers are not represented.
 */

typedef struct {
  size_t     depth;
  uint32_t   nCells;
  uint32_t   nAllocatedCells;
  cell2_t   *cells;
  uint32_t   nEdges;
  uint32_t   nAllocatedEdges;
  edge2_t   *edges;
  //cellInfo_t cellInfo;
} tree2_t;

void createTree2 (tree2_t *tree, size_t depth, uint32_t nCells, uint32_t nEdges) {
  tree->depth  = depth;
  tree->nCells = 0;
  tree->nEdges = 0;
  tree->nAllocatedCells = nCells;
  tree->nAllocatedEdges = nEdges;
  tree->cells = (cell2_t *) malloc(tree->nAllocatedCells * sizeof(cell2_t));
  tree->edges = (edge2_t *) malloc(tree->nAllocatedEdges * sizeof(edge2_t));
  //createCellInfo(&tree->cellInfo);
}

void freeTree2 (tree2_t *tree) {
  free(tree->cells);
  free(tree->edges);
  //freeCellInfo(&tree->cellInfo);
}

uint32_t addCell2 (tree2_t *tree, uint32_t edgeId, unsigned char nEdges) {
  assert(tree->nCells < tree->nAllocatedCells);
  assert(edgeId < tree->nEdges);
  tree->edges[edgeId].cellId = tree->nCells;
  tree->cells[tree->nCells].firstEdge = tree->nEdges;
  tree->cells[tree->nCells].nEdges    = nEdges;
  ++tree->nCells;
  tree->nEdges += nEdges;
  return tree->nCells-1;
}

uint32_t addEdge2 (tree2_t *tree, uint64_t cellId, unsigned char childId, unsigned char nucleotide, uint32_t sequence, uint32_t length) {
  cell2_t *cell = &tree->cells[cellId];
  assert(cellId < tree->nAllocatedCells);
  assert(cellId < tree->nCells);
  assert(childId < cell->nEdges);
  uint32_t edgeId = cell->firstEdge + childId;
  edge2_t *edge = &tree->edges[edgeId];
  edge->nucleotide = nucleotide;
  edge->sequence   = sequence;
  edge->length     = length;
  return edgeId;
}

/*
void setCellInfo (tree_t *tree, size_t cellId, size_t l, char *quality, unsigned int fileId) {
  ++tree->cells[cellId].counts[fileId];
  _setQuality(&tree->qualities, cellId, l, quality);
}
*/

bool isCellUnbranched (const tree2_t *tree, uint64_t cellId) {
  cell2_t *cell = &tree->cells[cellId];
  if (cell->nEdges == 0) {
    return true;
  }
  if (cell->nEdges > 1) {
    return false;
  }
  return isCellUnbranched(tree, tree->edges[cell->firstEdge].cellId);
}

#endif
