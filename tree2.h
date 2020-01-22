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
  uint32_t      cellId;
  char         *quality;
  count_t      *counts;
} cellInfo_t;

void createEmptyCellInfo (cellInfo_t *cellInfo) {
  cellInfo->cellId = NO_INFO;
}

bool isEmptyCellInfo (cellInfo_t *cellInfo) {
  return (cellInfo->cellId == NO_INFO);
}

void createCellInfo (cellInfo_t *cellInfo, size_t readSize, size_t nSamples) {
  ++readSize;
  cellInfo->quality = (char *)    malloc(readSize * sizeof(char));
  cellInfo->counts  = (count_t *) malloc(nSamples * sizeof(count_t));
}

void freeCellInfo (cellInfo_t *cellInfo) {
  free(cellInfo->quality);
  free(cellInfo->counts);
}

void setCellInfo (cellInfo_t *cellInfo, uint32_t cellId, char *quality, size_t readSize, count_t *counts, unsigned int nSamples) {
  cellInfo->cellId = cellId;
  ++readSize;
  memcpy(cellInfo->quality, quality, readSize * sizeof(char));
  memcpy(cellInfo->counts,  counts,  nSamples * sizeof(count_t));
}

void printCellInfo (cellInfo_t *cellInfo) {
  printf("cell id: %" PRIu32 ", qual: %s, counts: %p\n", cellInfo->cellId, cellInfo->quality, cellInfo->counts);
}

typedef struct {
  size_t nAllocatedCellInfos;
  size_t nCellInfos;
  cellInfo_t *cellInfos;
  cellInfo_t *currentCellInfo;
} cellInfos_t;

void createCellInfos (cellInfos_t *cellInfos, size_t nInfos, size_t readSize, size_t nSamples) {
  cellInfos->nAllocatedCellInfos = nInfos;
  cellInfos->nCellInfos = 0;
  cellInfos->cellInfos = (cellInfo_t *) malloc((nInfos+1) * sizeof(cellInfo_t));  
  for (unsigned int i = 0; i < nInfos; ++i) {
    createCellInfo(&cellInfos->cellInfos[i], readSize, nSamples);
  }
  createEmptyCellInfo(&cellInfos->cellInfos[nInfos]);
  cellInfos->currentCellInfo = cellInfos->cellInfos;
}

void freeCellInfos (cellInfos_t *cellInfos) {
  for (unsigned int i = 0; i < cellInfos->nCellInfos; ++i) {
    freeCellInfo(&cellInfos->cellInfos[i]);
  }
  free(cellInfos->cellInfos);
}

void addCellInfo (cellInfos_t *cellInfos, uint32_t cellId, char *quality, size_t readSize, count_t *counts, unsigned int nSamples) {
  assert(cellInfos->nCellInfos < cellInfos->nAllocatedCellInfos);
  setCellInfo(&cellInfos->cellInfos[cellInfos->nCellInfos], cellId, quality, readSize, counts, nSamples);
  ++cellInfos->nCellInfos;
}

cellInfo_t *getCellInfo (cellInfos_t *cellInfos, uint32_t cellId) {
  for (; (! isEmptyCellInfo(cellInfos->currentCellInfo)) && (cellInfos->currentCellInfo->cellId < cellId); ++cellInfos->currentCellInfo) {
  }
  return (cellInfos->currentCellInfo->cellId == cellId)? cellInfos->currentCellInfo: NULL;
}


typedef struct {
  uint32_t      firstEdge;
  unsigned char nEdges;
} cell2_t;

bool isTerminal (cell2_t *cell) {
  return (cell->nEdges == 0);
}

void printCell2(cell2_t *cell) {
  printf("(%p) %" PRIu32 "/%u", cell, cell->firstEdge, cell->nEdges);
}


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
  edge_t    *edges;
  cellInfos_t cellInfos;
} tree2_t;

void createTree2 (tree2_t *tree, size_t depth, uint32_t nCells, uint32_t nEdges, uint32_t nQualities) {
  printf("Creating new tree with depth %zu, %" PRIu32 " cells, %" PRIu32 " edges, %" PRIu32 " qualities\n", depth, nCells, nEdges, nQualities);
  tree->depth  = depth;
  tree->nCells = 0;
  tree->nEdges = 0;
  tree->nAllocatedCells = nCells;
  tree->nAllocatedEdges = nEdges;
  tree->cells = (cell2_t *) malloc(tree->nAllocatedCells * sizeof(cell2_t));
  tree->edges = (edge_t *) malloc(tree->nAllocatedEdges * sizeof(edge_t));
  createCellInfos(&tree->cellInfos, nQualities, tree->depth, parameters->nReadsFiles);
}

void freeTree2 (tree2_t *tree) {
  free(tree->cells);
  free(tree->edges);
  freeCellInfos(&tree->cellInfos);
}

edge_t *getFirstEdge2 (const tree2_t *tree, const cell2_t *cell) {
  return &tree->edges[cell->firstEdge];
}

uint32_t addCell2 (tree2_t *tree, uint32_t edgeId, unsigned char nEdges) {
  //printf("  Adding cell #%" PRIu32 "\n", tree->nCells);
  assert(tree->nCells < tree->nAllocatedCells);
  if (edgeId != NO_INFO) {
    assert(edgeId < tree->nEdges);
    tree->edges[edgeId].cellId = tree->nCells;
  }
  tree->cells[tree->nCells].firstEdge = tree->nEdges;
  tree->cells[tree->nCells].nEdges    = nEdges;
  ++tree->nCells;
  tree->nEdges += nEdges;
  return tree->nCells-1;
}

uint32_t addEdge2 (tree2_t *tree, uint32_t cellId, unsigned char childId, uint32_t sequence, uint32_t length) {
  assert(cellId < tree->nAllocatedCells);
  cell2_t *cell = &tree->cells[cellId];
  assert(childId < cell->nEdges);
  uint32_t edgeId = cell->firstEdge + childId;
  //printf("Adding edge2 with cell %" PRIu32 ", child %u, #cells %" PRIu32 "/%" PRIu32 ", #edges %" PRIu32 "/%" PRIu32 "\n", cellId, childId, tree->nCells, tree->nAllocatedCells, tree->nEdges, tree->nAllocatedEdges);
  //printCell2(cell); printf("\n"); fflush(stdout);
  assert(edgeId < tree->nAllocatedEdges);
  edge_t *edge    = &tree->edges[edgeId];
  edge->sequence  = sequence;
  edge->length    = length;
  return edgeId;
}

/*
void setCellInfo (tree_t *tree, size_t cellId, size_t l, char *quality, unsigned int fileId) {
  ++tree->cells[cellId].counts[fileId];
  _setQuality(&tree->qualities, cellId, l, quality);
}
*/

bool isCell2Unbranched (const tree2_t *tree, uint32_t cellId) {
  cell2_t *cell = &tree->cells[cellId];
  if (cell->nEdges == 0) {
    return true;
  }
  if (cell->nEdges > 1) {
    return false;
  }
  return isCell2Unbranched(tree, tree->edges[cell->firstEdge].cellId);
}

void _copyTree (tree2_t *tree2, uint32_t cellId2, const tree_t *tree, uint64_t cellId, size_t depth) {
  cell_t *cell = &tree->cells[cellId];
  edge_t *edge;
  char *quality;
  uint32_t edgeId2, newCellId2;
  unsigned short childId = 0;
  if ((quality = findQuality(&tree->qualities, cellId)) != NULL) {
    addCellInfo(&tree2->cellInfos, cellId2, quality, depth, cell->counts, parameters->nReadsFiles);
  }
  //printf("  Copying cell %" PRIu64 " to %" PRIu32 "/%" PRIu32 "\n", cellId, cellId2, tree2->nAllocatedCells);
  for (unsigned short nt = 0; nt < N_NUCLEOTIDES; ++nt) {
    edge = &cell->edges[nt]; 
    if (isSetEdge(edge)) {
      edgeId2 = addEdge2(tree2, cellId2, childId, edge->sequence, edge->length);
      newCellId2 = addCell2(tree2, edgeId2, getNChildren(&tree->cells[edge->cellId]));
      //printf("    edge #%i, edge id #%" PRIu32 "/%" PRIu32 ", cell id #%" PRIu32 "/%" PRIu32 "\n", nt, edgeId2, tree2->nAllocatedEdges, newCellId2, tree2->nAllocatedCells);
      _copyTree(tree2, newCellId2, tree, edge->cellId, depth + edge->length);
      ++childId;
    }
  }
}

void copyTree (tree2_t *tree2, const tree_t *tree) {
  createTree2(tree2, tree->depth, tree->nCells, tree->nEdges, tree->qualities.nQualities);
  for (uint32_t cellId = 0; cellId < N_TREE_BASE; ++cellId) {
    addCell2(tree2, NO_INFO, getNChildren(&tree->cells[cellId]));
  }
  tree2->nCells = N_TREE_BASE;
  for (uint32_t cellId = 0; cellId < N_TREE_BASE; ++cellId) {
    _copyTree(tree2, cellId, tree, cellId, TREE_BASE_SIZE);
  }
}

cellInfo_t *getCellInfoTree (tree2_t *tree, uint32_t cellId) {
  return getCellInfo(&tree->cellInfos, cellId);
}


#endif
