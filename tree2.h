#ifndef TREE2_H
#define TREE2_H

#include <assert.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "constants.h"
#include "parameters.h"
#include "quality.h"
#include "edge.h"
#include "cell.h"
#include "tree.h"

#define NO_INFO ((uint32_t) -1)

/**
 * CellInfo stores informations on compressed tree cells:
 *  - the id of the cell
 *  - the quality
 *  - the counts
 */
typedef struct {
  uint32_t      cellId;
  char         *quality;
  count_t      *counts;
} cellInfo_t;

void createEmptyCellInfo (cellInfo_t *cellInfo) {
  cellInfo->cellId = NO_INFO;
}

bool isEmptyCellInfo (const cellInfo_t *cellInfo) {
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

/**
 * CellInfos store all the cell infos in a vector.
 * The infos a supposed to be read sequentially, from the first cell to the last (possibly skipping some)...
 */
typedef struct {
  size_t nAllocatedCellInfos;
  size_t nCellInfos;
  cellInfo_t *cellInfos;
} cellInfos_t;

typedef struct {
  cellInfo_t *currentCellInfo;
} cellVisitor_t;

void createCellVisotor (cellVisitor_t *cellVisitor, cellInfos_t *cellInfos) {
  cellVisitor->currentCellInfo = cellInfos->cellInfos;
}

void createCellInfos (cellInfos_t *cellInfos, size_t nInfos, size_t readSize, size_t nSamples) {
  cellInfos->nAllocatedCellInfos = nInfos;
  cellInfos->nCellInfos = 0;
  cellInfos->cellInfos = (cellInfo_t *) malloc((nInfos+1) * sizeof(cellInfo_t));  
  for (unsigned int i = 0; i < nInfos; ++i) {
    createCellInfo(&cellInfos->cellInfos[i], readSize, nSamples);
  }
  createEmptyCellInfo(&cellInfos->cellInfos[nInfos]);
}

void freeCellInfos (cellInfos_t *cellInfos) {
  for (unsigned int i = 0; i < cellInfos->nAllocatedCellInfos; ++i) {
    freeCellInfo(&cellInfos->cellInfos[i]);
  }
  free(cellInfos->cellInfos);
}

void addCellInfo (cellInfos_t *cellInfos, uint32_t cellId, char *quality, size_t readSize, count_t *counts, unsigned int nSamples) {
  assert(cellInfos->nCellInfos < cellInfos->nAllocatedCellInfos);
  //printf("Setting cell info #%" PRIu32 "\n", cellId);
  setCellInfo(&cellInfos->cellInfos[cellInfos->nCellInfos], cellId, quality, readSize, counts, nSamples);
  ++cellInfos->nCellInfos;
}

cellInfo_t *getCellInfo (uint32_t cellId, cellVisitor_t *cellVisitor) {
  //printf("Looking for cell info %" PRIu32 "\n", cellId); fflush(stdout);
  //printf("Current state is %p %" PRIu32 "\n", cellInfos->currentCellInfo, cellInfos->currentCellInfo->cellId); fflush(stdout);
  //printf("Current state is %p\n", cellInfos->cellInfos + cellInfos->nCellInfos); fflush(stdout);
  for (; (! isEmptyCellInfo(cellVisitor->currentCellInfo)) && (cellVisitor->currentCellInfo->cellId < cellId); ++cellVisitor->currentCellInfo) {
    //printf("Got %" PRIu32 "\n", cellInfos->currentCellInfo->cellId); fflush(stdout);
  }
  return (cellVisitor->currentCellInfo->cellId == cellId)? cellVisitor->currentCellInfo: NULL;
}

/**
 * A cell is:
 * - the index of the first (child) edge in the tree's vector of edges
 * - the number of edges
 */
typedef struct {
  uint32_t      firstEdge;
  unsigned char nEdges;
} cell2_t;

unsigned char getNChildren2 (const cell2_t *cell) {
  return cell->nEdges;
}

bool isTerminal (const cell2_t *cell) {
  return (cell->nEdges == 0);
}

void printCell2 (const cell2_t *cell) {
  printf("(%p) %" PRIu32 "/%u", cell, cell->firstEdge, cell->nEdges);
}


/******* Compressed tree type *******/
/**
 * A tree is stores all the prefixes of the reads.
 * It is:
 *   - a vector of edges
 *   - a vector of cells
 *   - a depth
 *   - the number of cells currently used
 *   - the number of allocated cells
 *   - the number of edges currently used
 *   - the number of allocated edges
 *   - a pointer to the vector of cellInfos
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
  //printf("  Adding cell #%" PRIu32 "/%" PRIu32 " after edge %" PRIu32 "/%" PRIu32 " (NO_INFO = %" PRIu32 ")\n", tree->nCells, tree->nAllocatedCells, edgeId, tree->nEdges, NO_INFO); fflush(stdout);
  assert(tree->nCells < tree->nAllocatedCells);
  if (edgeId != NO_INFO) {
    assert(edgeId < tree->nEdges);
    tree->edges[edgeId].cellId = tree->nCells;
  }
  tree->cells[tree->nCells].firstEdge = tree->nEdges;
  tree->cells[tree->nCells].nEdges    = nEdges;
  ++tree->nCells;
  assert(tree->nCells <= tree->nAllocatedCells);
  tree->nEdges += nEdges;
  assert(tree->nEdges <= tree->nAllocatedEdges);
  return tree->nCells-1;
}

uint32_t addEdge2 (tree2_t *tree, uint32_t cellId, unsigned char childId, uint32_t sequence, uint32_t length) {
  assert(cellId < tree->nAllocatedCells);
  cell2_t *cell = &tree->cells[cellId];
  assert(childId < cell->nEdges);
  uint32_t edgeId = cell->firstEdge + childId;
  //printf("Adding edge2 with cell %" PRIu32 ", child %u, #cells %" PRIu32 "/%" PRIu32 ", #edges %" PRIu32 "/%" PRIu32 "\n", cellId, childId, tree->nCells, tree->nAllocatedCells, tree->nEdges, tree->nAllocatedEdges);
  //printCell2(cell); printf("\n"); fflush(stdout);
  assert(edgeId < tree->nEdges);
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

/**
 * Copy a non-compressed tree to a compressed tree
 */
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

cellInfo_t *getCellInfoTree (uint32_t cellId, cellVisitor_t *cellVisitor) {
  return getCellInfo(cellId, cellVisitor);
}


#endif
