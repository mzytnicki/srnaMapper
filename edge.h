#ifndef EDGE_H
#define EDGE_H

/******* Edge type *******/
/**
 * An edge is the link from one cell to one of its children.
 * It also contains the sequence for a cell to the child.
 * The sequence is stored in a int, with the first nucleotides being the sequence, and the last ones the size
 * It is:
 *   - the sequence: mix of sequence and length
 *   - the id of the next cell
 */
typedef uint32_t sequence_t;

typedef struct {
  sequence_t sequence: EDGE_SEQ_LENGTH;
  sequence_t length:   EDGE_LENGTH_LENGTH;
  uint32_t cellId;
} edge_t;

void createEdge (edge_t *edge) {
  edge->length   = 0;
  edge->sequence = 0;
  edge->cellId   = NO_DATA;
}

bool isSetEdge (edge_t *edge) {
  return (edge->length > 0);
}

bool isFullEdge (edge_t *edge) {
  return (edge->length == MAX_EDGE_LENGTH);
}

void setEdge (edge_t *edge, sequence_t sequence, sequence_t length, uint64_t cellId) {
  edge->sequence = sequence;
  edge->length   = length;
  edge->cellId   = cellId;
}

void unsetEdge(edge_t *edge) {
  edge->sequence = 0;
  edge->length   = 0;
  edge->cellId   = NO_DATA;
}

void setCellIdEdge (edge_t *edge, uint64_t cellId) {
  edge->cellId = cellId;
}

void addEdgeNucleotide (edge_t *edge, unsigned short nucleotide) {
  assert(! isFullEdge(edge));
  //printf("    adding nucl. %u to %i @ %i\n", nucleotide, edge->sequence, edge->length);
  edge->sequence |= nucleotide << (NUCLEOTIDES_BITS * edge->length);
  //printf("      now: %i\n", edge->sequence);
  ++edge->length;
}

unsigned short getEdgeNucleotide (edge_t *edge, sequence_t length) {
  assert(length <= edge->length);
  assert(length <= MAX_EDGE_LENGTH);
  return getNucleotide(edge->sequence, length);
}

void edgeRemoveFirst (edge_t *edge, size_t length) {
  assert(length <= edge->length);
  edge->length -= length;
  edge->sequence >>= (length * NUCLEOTIDES_BITS);
}

void edgeRemoveLast (edge_t *edge, size_t newSize) {
  assert(newSize <= edge->length);
  edge->length = newSize;
  edge->sequence &= ((1 << (newSize * NUCLEOTIDES_BITS)) - 1);
}

void edgeSetCellId (edge_t *edge, uint64_t cellId) {
  edge->cellId = cellId;
}

void printEdge (edge_t *edge) {
  if (! isSetEdge(edge)) {
    printf("-X->");
    return;
  }
  printf("-");
  printSequence(edge->sequence, edge->length);
  printf("-> %" PRIu32, edge->cellId);
}

#endif
