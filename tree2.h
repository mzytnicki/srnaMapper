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

void printCellInfo (const cellInfo_t *cellInfo) {
  printf("cell id: %" PRIu32 ", qual: %s, counts (%p): ", cellInfo->cellId, cellInfo->quality, cellInfo->counts);
  for (unsigned int countId = 0; countId < parameters->nReadsFiles; ++countId) {
    printf("%" PRIu32 " ", cellInfo->counts[countId]);
  }
  printf("\n");
}

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
} cellInfos_t;

typedef size_t cellVisitor_t;

void clearCellVisitor (cellVisitor_t *cellVisitor) {
  *cellVisitor = 0;
}

void createCellInfos (cellInfos_t *cellInfos, size_t nInfos, size_t readSize, size_t nSamples) {
  cellInfos->nAllocatedCellInfos = nInfos;
  cellInfos->nCellInfos          = 0;
  cellInfos->cellIds             = (uint32_t *) malloc((nInfos+1) * sizeof(uint32_t));  
  cellInfos->counts              = (count_t *) malloc(nInfos * (parameters->nReadsFiles) * sizeof(count_t));  
  cellInfos->qualities           = (char **) malloc(nInfos * sizeof(char *));  
  cellInfos->cellIds[nInfos]     = NO_INFO;
}

void freeCellInfos (cellInfos_t *cellInfos) {
  for (size_t cellInfosId = 0; cellInfosId < cellInfos->nCellInfos; ++cellInfosId) {
    free(cellInfos->qualities[cellInfosId]);
  }
  free(cellInfos->cellIds);
  free(cellInfos->qualities);
  free(cellInfos->counts);
}

void addCellInfo (cellInfos_t *cellInfos, uint32_t cellId, char *quality, size_t readSize, count_t *counts, unsigned int nSamples) {
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
  ++cellInfos->nCellInfos;
}

//cellInfo_t *getCellInfo (uint32_t cellId, cellVisitor_t *cellVisitor) {
bool getCellInfo (const cellInfos_t *cellInfos, uint32_t cellId, cellVisitor_t *cellVisitor) {
  //printf("Looking for cell info %" PRIu32 "\n", cellId); fflush(stdout);
  //printf("Current state is %zu -> %" PRIu32 "\n", cellVisitor, cellInfos->cellIds[*cellVisitor]); fflush(stdout);
  //printf("Current state is %p\n", cellInfos->cellInfos + cellInfos->nCellInfos); fflush(stdout);
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


/*
void printCellInfos (const cellInfos_t *cellInfos) {
  for (size_t cellInfoId = 0; cellInfoId < cellInfos->nCellInfos; ++cellInfoId) {
    printCellInfo(&cellInfos->cellInfos[cellInfoId]);
  }
}
*/

#define EDGE2_SEQ_START_LENGTH  26
#define MAX_EDGE2_SEQ_START     (1 << EDGE2_SEQ_START_LENGTH)
#define EDGE2_LENGTH_LENGTH      6
#define MAX_EDGE2_LENGTH        (1 << EDGE2_LENGTH_LENGTH)
#define EDGE_SEQ_DIRECT_DEPTH   (EDGE2_SEQ_START_LENGTH / NUCLEOTIDES_BITS)

typedef struct {
  uint32_t sequenceStart: EDGE_SEQ_LENGTH;
  uint32_t length:        EDGE_LENGTH_LENGTH;
  uint32_t cellId;
} edge2_t;

void createEdge2(edge2_t *edge) {
  edge->length = 0;
}

bool isSetEdge2(const edge2_t *edge) {
  return (edge->length == 0);
}

void printEdge2 (edge2_t *edge) {
  if (! isSetEdge2(edge)) {
    printf("-X->");
    return;
  }
  printf(" (%u)-> %" PRIu32, edge->length, edge->cellId);
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

#define EDGE_SEQ_ALLOC 0x10000

typedef struct {
  size_t     depth;
  uint32_t   nCells;
  uint32_t   nAllocatedCells;
  cell2_t   *cells;
  uint32_t   nEdges;
  uint32_t   nAllocatedEdges;
  edge_t    *edges;
  /*
  edge2_t   *edges;
  uint32_t   edgeSequenceAllocated;
  uint32_t   edgeSequenceLength;
  char      *edgeSequence;
  */
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
  /*
  tree->edges = (edge2_t *) malloc(tree->nAllocatedEdges * sizeof(edge2_t));
  tree->edgeSequenceLength    = 0;
  tree->edgeSequenceAllocated = EDGE_SEQ_ALLOC;
  tree->edgeSequence = (char *) calloc(tree->edgeSequenceAllocated, sizeof(char));
  */
  createCellInfos(&tree->cellInfos, nQualities, tree->depth, parameters->nReadsFiles);
}

void freeTree2 (tree2_t *tree) {
  free(tree->cells);
  free(tree->edges);
  //free(tree->edgeSequence);
  freeCellInfos(&tree->cellInfos);
}

//edge2_t *getFirstEdge2 (const tree2_t *tree, const cell2_t *cell) {
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

/*
#define N_NUCLEOTIDES_PER_BYTE 4
#define EDGE_SEQUENCE_SHIFT 2
#define EDGE_SEQUENCE_MASK  3

uint32_t getEdgeSequenceArrayOffset (uint32_t position) {
  return position >> EDGE_SEQUENCE_SHIFT;
}

uint32_t getEdgeSequenceCharOffset (uint32_t position) {
  return (N_NUCLEOTIDES_PER_BYTE - (position & EDGE_SEQUENCE_MASK) - 1) * NUCLEOTIDES_BITS;
}

int getEdge2Nucleotide(const tree2_t *tree, const edge2_t *edge, uint32_t position) {
  assert(position <= edge->length);
  if (edge->length < EDGE_SEQ_DIRECT_DEPTH) {
    return getNucleotide(edge->sequenceStart, position);
  }
  return (tree->edgeSequence[getEdgeSequenceArrayOffset(edge->sequenceStart + position)] >> getEdgeSequenceCharOffset(edge->sequenceStart + position)) & NUCLEOTIDE_MASK;
}

void printEdge2Sequence (const tree2_t *tree, edge2_t *edge) {
  if (edge->length < EDGE_SEQ_DIRECT_DEPTH) {
    printSequence(edge->sequenceStart, edge->length);
    return;
  }
  for (uint32_t position = 0; position < edge->length; ++position) {
    printf("%c", "ACGT"[getEdge2Nucleotide(tree, edge, position)]);
  }
}

void setEdgeSequence(tree2_t *tree, uint32_t position, unsigned char nucleotide) {
  tree->edgeSequence[getEdgeSequenceArrayOffset(position)] |= (nucleotide << getEdgeSequenceCharOffset(position));
}

uint32_t getEdgeSequencesSize (uint32_t nNucleotides) {
  return (nNucleotides / N_NUCLEOTIDES_PER_BYTE);
}

void copyEdgeSequence (tree2_t *tree, uint32_t sequence, uint32_t length) {
  assert(length >= EDGE_SEQ_DIRECT_DEPTH);
  if (getEdgeSequencesSize(tree->edgeSequenceLength + length) >= tree->edgeSequenceAllocated) {
    tree->edgeSequenceAllocated *= 2;
    tree->edgeSequence = (char *) realloc(tree->edgeSequence, tree->edgeSequenceAllocated * sizeof(char));
    if (tree->edgeSequence == NULL) {
      fprintf(stderr, "Error!  Cannot allocate a vector of size %" PRIu32 " for the edge sequences.\nExiting.\n", tree->edgeSequenceAllocated);
    }
    printf("Allocating a vector of size %" PRIu32 " for the edge sequences.\n", tree->edgeSequenceAllocated); fflush(stdout);
    if (tree->edgeSequenceAllocated >= MAX_EDGE2_SEQ_START) {
      fprintf(stderr, "Error!  Trying to allocate an edge sequence of size %" PRIu32 " > %u.\nExiting.\n", tree->edgeSequenceAllocated, MAX_EDGE2_SEQ_START);
      exit(EXIT_FAILURE);
    }
  }
  for (uint32_t sequenceId = 0; sequenceId < length; ++sequenceId) {
    char c = sequence & NUCLEOTIDE_MASK;
    setEdgeSequence(tree, tree->edgeSequenceLength, c);
    ++tree->edgeSequenceLength;
    sequence >>= NUCLEOTIDES_BITS;
  }
}
*/

uint32_t addEdge2 (tree2_t *tree, uint32_t cellId, unsigned char childId, uint32_t sequence, uint32_t length) {
  assert(cellId < tree->nAllocatedCells);
  cell2_t *cell = &tree->cells[cellId];
  assert(childId < cell->nEdges);
  uint32_t edgeId = cell->firstEdge + childId;
  //printf("Adding edge2 with cell %" PRIu32 ", child %u, #cells %" PRIu32 "/%" PRIu32 ", #edges %" PRIu32 "/%" PRIu32 "\n", cellId, childId, tree->nCells, tree->nAllocatedCells, tree->nEdges, tree->nAllocatedEdges);
  //printCell2(cell); printf("\n"); fflush(stdout);
  assert(edgeId < tree->nEdges);
  assert(edgeId < tree->nAllocatedEdges);
  //edge2_t *edge       = &tree->edges[edgeId];
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

bool getCellInfoTree (const tree2_t *tree, uint32_t cellId, cellVisitor_t *cellVisitor) {
  return getCellInfo(&tree->cellInfos, cellId, cellVisitor);
}

/**
 * Print the last nucleotides (close to the leaves) of the tree.
 */
void __printTree2 (const tree2_t *tree, char *read, size_t readPos, uint32_t cellId, cellVisitor_t *cellVisitor) {
  edge_t *edge;
  //edge2_t *edge;
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
  char *read = (char *) malloc((tree->depth+1) * sizeof(char));
  cellVisitor_t cellVisitor;
  clearCellVisitor(&cellVisitor);
  read[tree->depth] = 0;
  _printTree2(tree, read, 0, 0, &cellVisitor);
  free(read);
  //printf("Printing cell info\n");
  //printCellInfos(&tree->cellInfos);
}


#endif
