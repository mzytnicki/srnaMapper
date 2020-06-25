#ifndef CELL2_H
#define CELL2_H

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

#endif
