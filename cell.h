#ifndef CELL_H
#define CELL_H

/******* Cell type *******/
/**
 * A cell is a prefix of the reads tree.
 * It is:
 *   - 4 children: an array of 4 indices of the children array in the tree structure
 *   - the read counts: an array of int
 */

typedef struct {
  edge_t edges [N_NUCLEOTIDES];
  count_t *counts;
} cell_t;

void createCell (cell_t *cell) {
  for (unsigned short i = 0; i < N_NUCLEOTIDES; ++i) {
    createEdge(&cell->edges[i]);
  }
  //printf("%u %zu %p\n", parameters->nReadsFiles, sizeof(count_t), cell);
  cell->counts = (count_t *) calloc(parameters->nReadsFiles, sizeof(count_t));
}

void freeCell (cell_t *cell) {
  free(cell->counts);
}

void addEdge (cell_t *cell, sequence_t sequence, sequence_t length, uint64_t childId) {
  unsigned short nucleotide = sequence & NUCLEOTIDE_MASK;
  setEdge(&cell->edges[nucleotide], sequence, length, childId);
}

edge_t *getFirstEdge (cell_t *cell) {
  edge_t *edge;
  for (unsigned short edgeId = 0; edgeId < N_NUCLEOTIDES; ++edgeId) {
    edge = &cell->edges[edgeId];
    if (isSetEdge(edge)) {
      return edge;
    }
  }
  return NULL;
}

/**
 * Split an edge into two.  Add the remaining part of the edge to a new cell.
 * When newEdgeId == N_NUCLEOTIDES, the nucleotide at the split point is not known
 */
void splitEdge (edge_t *edge, cell_t *newCell, uint64_t newCellId, sequence_t length, unsigned short newEdgeId) {
  edge_t *newEdge;
  //printf("   Splitting edge %p %i into %u -> %" PRIu64 " at point %zu and id %i\n", edge, edge->sequence, edge->length, edge->cellId, length, newEdgeId);
  assert(edge->cellId != NO_DATA);
  assert(length < edge->length);
  assert(newCellId != NO_DATA);
  if (newEdgeId == N_NUCLEOTIDES) newEdgeId = getEdgeNucleotide(edge, length);
  newEdge = &newCell->edges[newEdgeId];
  newEdge->sequence = edge->sequence;
  newEdge->length = edge->length;
  edgeRemoveFirst(newEdge, length);
  edgeRemoveLast(edge, length);
  newEdge->cellId = edge->cellId;
  edge->cellId = newCellId;
  //printf("   split edge %p into %i (%u) -> %" PRIu64 " and %p %i (%u) -> %" PRIu64 "\n", edge, edge->sequence, edge->length, edge->cellId, newEdge, newEdge->sequence, newEdge->length, newEdge->cellId);
}

unsigned short getNChildren (cell_t *cell) {
  unsigned short nChildren = 0;
  for (unsigned int edgeId = 0; edgeId < N_NUCLEOTIDES; ++edgeId) {
    if (isSetEdge(&cell->edges[edgeId])) {
      ++nChildren;
    }
  }
  return nChildren;
}

/*
bool isTerminal (cell_t *cell) {
  for (unsigned int edgeId = 0; edgeId < N_NUCLEOTIDES; ++edgeId) {
    if (isSetEdge(&cell->edges[edgeId])) {
      return false;
    }
  }
  return true;
}
*/

void printCell (cell_t *cell) {
  for (unsigned short i = 0; i < N_NUCLEOTIDES; ++i) {
    if (isSetEdge(&cell->edges[i])) {
      printEdge(&cell->edges[i]);
      printf("  ");
    }
  }
}

#endif
