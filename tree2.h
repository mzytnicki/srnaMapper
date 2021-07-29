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
#include "cellInfo.h"
#include "cell2.h"

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

#define EDGE_SEQ_ALLOC 0x10000

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
  tree->cells = (cell2_t *) mallocOrDie(tree->nAllocatedCells * sizeof(cell2_t));
  tree->edges = (edge_t *) mallocOrDie(tree->nAllocatedEdges * sizeof(edge_t));
  createCellInfos(&tree->cellInfos, nQualities);
}

void freeTree2 (tree2_t *tree) {
  free(tree->cells);
  free(tree->edges);
  //free(tree->edgeSequence);
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
  edge_t *edge       = &tree->edges[edgeId];
  edge->length        = length;
  edge->sequence = sequence;
  /*
  if (length >= EDGE_SEQ_DIRECT_DEPTH) {
    edge->sequenceStart = tree->edgeSequenceLength;
    copyEdgeSequence(tree, sequence, length);
  }
  else {
    edge->sequenceStart = sequence;
  }
  */
  //printf("Compare sequence (size %" PRIu32 "): ", length);
  //printSequence(sequence, length);
  //printf(" with edge ");
  //printEdge2Sequence(tree, edge); printf("\n"); fflush(stdout);
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
    //printf("  Adding cell info @ %"PRIu32"\n", cellId2);
    addCellInfo(&tree2->cellInfos, cellId2, quality, depth, cell->counts);
  }
  //if (depth > TREE_BASE_SIZE) printf("  Copying cell %" PRIu64 " to %" PRIu32 "/%" PRIu32 "\n", cellId, cellId2, tree2->nAllocatedCells);
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

bool getCellInfoTree (const tree2_t *tree, uint32_t cellId, cellVisitor_t *cellVisitor) {
  return getCellInfo(&tree->cellInfos, cellId, cellVisitor);
}

/**
 * Print the last nucleotides (close to the leaves) of the tree.
 */
void __printTree2 (const tree2_t *tree, char *read, size_t readPos, uint32_t cellId, cellVisitor_t *cellVisitor) {
  edge_t *edge;
  cell2_t *cell = &tree->cells[cellId];
  //printf("Read: %s at %"PRIu32 " with %" PRIu32 " children at pos %zu/%zu\n", read+tree->depth-readPos, cellId, cell->nEdges, readPos, tree->depth); fflush(stdout);
  if (getCellInfoTree(tree, cellId, cellVisitor)) {
    printf("%s\n", read+tree->depth-readPos);
    //assert(strlen(read+tree->depth-readPos) == strlen(cell->quality));
  }
  for (uint32_t edgeId = 0; edgeId < cell->nEdges; ++edgeId) {
    edge = &tree->edges[cell->firstEdge + edgeId];
    //printf("Following to node %" PRIu32 " -> %" PRIu32 " and length %" PRIu32 "\n", edgeId, edge->cellId, edge->length); fflush(stdout);
    for (size_t i = 0; i < edge->length; ++i) {
      assert(1 + readPos + i <= tree->depth);
      //read[tree->depth-readPos-1-i] = DNA5_TO_CHAR[getEdge2Nucleotide(tree, edge, i)];
      read[tree->depth-readPos-1-i] = DNA5_TO_CHAR[getNucleotide(edge->sequence, i)];
    }
    __printTree2(tree, read, readPos + edge->length, edge->cellId, cellVisitor);
  }
}

/**
 * Print the first nucleotides (close to the root) of the tree, then call __printTree.
 */
void _printTree2 (const tree2_t *tree, char *read, size_t readPos, uint32_t cellId, cellVisitor_t *cellVisitor) {
  if (readPos == TREE_BASE_SIZE) {
    __printTree2(tree, read, readPos, cellId, cellVisitor);
    return;
  }
  assert(cellId < N_TREE_BASE);
  cellId <<= NUCLEOTIDES_BITS;
  for (unsigned short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    read[tree->depth-readPos-1] = DNA5_TO_CHAR[nucleotide];
    _printTree2(tree, read, readPos+1, cellId+nucleotide, cellVisitor);
  }
}

/**
 * Open/close file, allocate the memory, and call _printTree
 */
void printTree2 (const tree2_t *tree) {
  printf("Printing tree\n");
  char *read = (char *) mallocOrDie((tree->depth+1) * sizeof(char));
  cellVisitor_t cellVisitor;
  clearCellVisitor(&cellVisitor);
  read[tree->depth] = 0;
  _printTree2(tree, read, 0, 0, &cellVisitor);
  free(read);
  //printf("Printing cell info\n");
  //printCellInfos(&tree->cellInfos);
}


#endif
