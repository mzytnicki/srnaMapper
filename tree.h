#ifndef TREE_H
#define TREE_H

#include "constants.h"
#include "parameters.h"
#include "helper.h"
#include "edge.h"
#include "cell.h"
#include "quality.h"

/******* Tree type *******/
/**
 * A tree is stores all the prefixes of the reads.
 * It is:
 *   - a depth
 *   - an array of cells
 *   - the number of cells
 *   - the size of the array
 *   - an instance of the quality structure
 * The first N_TREE_BASE cells correspond to the all the combinations of TREE_BASE_SIZE-mers.
 * In other words, the 1-mers, 2-mers, ..., TREE_BASE_SIZE-1-mers are not represented.
 */

typedef struct {
  size_t    depth;
  cell_t   *cells;
  uint64_t  nCells;
  uint64_t  nEdges;
  uint64_t  nAllocated;
  quality_t qualities;
} tree_t;

void createTree (tree_t *tree) {
  tree->nAllocated = INIT_N_CELLS;
  tree->cells      = (cell_t *) malloc(tree->nAllocated * sizeof(cell_t));
  tree->depth      = TREE_BASE_SIZE;
  tree->nCells     = N_TREE_BASE;
  tree->nEdges     = 0;
  for (size_t i = 0; i < N_TREE_BASE; ++i) {
    createCell(&tree->cells[i]);
  }
  createQualities(&tree->qualities);
}

void _computeTreeStats (const tree_t *tree, unsigned int **stats, unsigned int *statsSum, unsigned int *branchSizes, unsigned int branchSize, cell_t *cell, size_t depth, unsigned int *nNodes, unsigned int *nQualities) {
  edge_t *edge;
  size_t length;
  unsigned int c = 0;
  for (short edgeId = 0; edgeId < N_NUCLEOTIDES; ++edgeId) {
    edge = &cell->edges[edgeId];
    if (edge->cellId != NO_DATA) {
      length = edge->length;
      ++c;
      ++(*nNodes);
      if (findQualityId(&tree->qualities, edge->cellId) != NO_QUALITY) {
        ++(*nQualities);
      }
    }
  }
  if (c == 1) {
    branchSize += length;
  }
  else {
    ++branchSizes[branchSize];
    branchSize = 0;
  }
  for (short edgeId = 0; edgeId < N_NUCLEOTIDES; ++edgeId) {
    edge = &cell->edges[edgeId];
    if (edge->cellId != NO_DATA) {
      _computeTreeStats(tree, stats, statsSum, branchSizes, branchSize, &tree->cells[edge->cellId], depth+edge->length, nNodes, nQualities);
    }
  }
  //printf("At depth %zu/%zu with %i nucleotides\n", depth, tree->depth, c); fflush(stdout);
  ++stats[depth][c];
  ++statsSum[c];
}

void computeTreeStats (const tree_t *tree) {
  //unsigned int **stats = (unsigned int **) malloc((tree->depth+1) * N_NUCLEOTIDES * sizeof(unsigned int *));
  unsigned int **stats = (unsigned int **) malloc((tree->depth+1) * sizeof(unsigned int *));
  unsigned int statsSum[N_NUCLEOTIDES] = { 0, 0, 0, 0 };
  unsigned int *branchSizes = (unsigned int *) calloc(tree->depth+1, sizeof(unsigned int));
  unsigned int nNodes = 0;
  unsigned int nQualities = 0;
  unsigned int s;
  for (size_t depth = 0; depth <= tree->depth; ++depth) {
    stats[depth] = (unsigned int *) calloc(N_NUCLEOTIDES + 1, sizeof(unsigned int));
  }
  for (size_t cellId = 0; cellId < N_TREE_BASE; ++cellId) {
    _computeTreeStats(tree, stats, statsSum, branchSizes, 0, &tree->cells[cellId], TREE_BASE_SIZE, &nNodes, &nQualities);
  }
  puts("Stats on tree");
  for (size_t depth = TREE_BASE_SIZE; depth <= tree->depth; ++depth) {
    printf("%zu:", depth);
    s = 0;
    for (size_t nChildren = 0; nChildren < N_NUCLEOTIDES; ++nChildren) {
      s += stats[depth][nChildren];
    }
    if (s != 0) {
      for (size_t nChildren = 0; nChildren < N_NUCLEOTIDES; ++nChildren) {
        printf("\t%zu:%u (%.f%%)", nChildren, stats[depth][nChildren], ((float) stats[depth][nChildren]) / s * 100);
      }
    }
    printf("\n");
  }
  printf("sum:");
  s = 0;
  for (size_t nChildren = 0; nChildren < N_NUCLEOTIDES; ++nChildren) {
    s += statsSum[nChildren];
  }
  for (size_t nChildren = 0; nChildren < N_NUCLEOTIDES; ++nChildren) {
    printf("\t%zu:%u (%.f%%)", nChildren, statsSum[nChildren], ((float) statsSum[nChildren]) / s * 100);
  }
  printf("\n");
  puts("Stats on branch sizes");
  for (size_t depth = 0; depth <= tree->depth; ++depth) {
    printf("%zu: %u\n", depth, branchSizes[depth]);
  }
  printf("%u nodes with %u qualities (%.f%%)\n", nNodes, nQualities, ((float) nQualities) / nNodes * 100);
  for (size_t depth = 0; depth <= tree->depth; ++depth) {
    free(stats[depth]);
  }
  free(stats);
  free(branchSizes);
}

void freeTree (tree_t *tree) {
  for (uint64_t i = 0; i < tree->nCells; ++i) {
    freeCell(&tree->cells[i]);
  }
  free(tree->cells);
  freeQualities(&tree->qualities);
  tree->nAllocated = 0;
  tree->depth = 0;
  tree->nCells = 0;
}

void setQuality (tree_t *tree, size_t cellId, size_t l, char *quality, unsigned int fileId) {
  ++tree->cells[cellId].counts[fileId];
  _setQuality(&tree->qualities, cellId, l, quality);
}

uint64_t addCell (tree_t *tree) {
  if (tree->nCells == tree->nAllocated) {
    tree->nAllocated *= 2;
    //printf("reallocating cells: %zu to %zu...\n", tree->nCells, tree->nAllocated);
    if ((tree->cells = (cell_t *) realloc(tree->cells, tree->nAllocated * sizeof(cell_t))) == NULL) {
      printf("Cannot allocate memory for tree of size %" PRIu64 ".\nExiting.\n", tree->nAllocated);
      exit(EXIT_FAILURE);
    }
  }
  //printf("ncells: %zu / %zu\n", tree->nCells, tree->nAllocated);
  createCell(&tree->cells[tree->nCells]);
  ++tree->nCells;
  return tree->nCells-1;
}

/**
 * Create a new node, split an edge into two, and add the remaining part of the edge to the new cell.
 * The current edge keeps the first part of the sequence.
 */
uint64_t splitEdgeTree (tree_t *tree, uint64_t cellId, unsigned short edgeId, sequence_t length, unsigned short newEdgeId) {
  // Adding cell can trigger reallocation and invalidate edge pointers
  uint64_t newCellId = addCell(tree);
  splitEdge(&tree->cells[cellId].edges[edgeId], &tree->cells[newCellId], newCellId, length, newEdgeId);
  return newCellId;
}

bool isCellUnbranched (const tree_t *tree, cell_t *cell) {
  bool foundOneChild = false;
  unsigned short child = 0;
  //printf("\t\tCell unbranched: ");
  //printCell(cell);
  //printf("\n");
  for (unsigned int edgeId = 0; edgeId < N_NUCLEOTIDES; ++edgeId) {
    if (isSetEdge(&cell->edges[edgeId])) {
      //printEdge(&cell->edges[edgeId]);
      if (foundOneChild) {
        //printf("\t\t\tNo\n");
        return false;
      }
      foundOneChild = true;
      child = edgeId;
    }
  }
  if (foundOneChild) {
    return isCellUnbranched(tree, &tree->cells[cell->edges[child].cellId]);
  }
  //printf("\t\t\tYes\n");
  return true;
}

/**
 * Add the rest of the sequence to the tree and extend the edge.
 * Explore the tree by following the sequence.
 * Do not add anything the tree ("addSequenceAdd" does it).
 * Start at node cellId and character sequenceId.
 */
uint64_t addSequenceAdd (tree_t *tree, uint64_t cellId, char *binarySequence, int sequenceId, edge_t *edge) {
  unsigned short nucleotide;
  assert(sequenceId >= 0);
  assert(cellId < tree->nCells);
  //printf("Adding sequence %s @%" PRIu64 "\n", sequence, cellId);
  for (; sequenceId >= 0; --sequenceId) {
    nucleotide = binarySequence[sequenceId];
    if (edge == NULL) {
      ++tree->nEdges;
      edge = &tree->cells[cellId].edges[nucleotide];
    }
    else if (edge->length == 0) {
      ++tree->nEdges;
    }
    //printf("  seq id: %d, edge: %p (%i), nucleotide: %c\n", sequenceId, edge, edge->length, DNA5_TO_CHAR[nucleotide]);
    addEdgeNucleotide(edge, nucleotide);
    if (edge->length == MAX_EDGE_LENGTH) {
      setCellIdEdge(edge, tree->nCells);
      // Adding cell may trigger reallocate and invalidate edge pointer
      cellId = addCell(tree);
      //printf("  to the end: %" PRIu64 "\n", cellId);
      //setCellIdEdge(edge, cellId);
      edge = NULL;
    }
  }
  // set node
  if (edge != NULL) {
      // Adding cell may trigger reallocate and invalidate edge pointer
      setCellIdEdge(edge, tree->nCells);
      cellId = addCell(tree);
  }
  //printf("  over with %" PRIu64 "\n", cellId);
  return cellId;
}

/**
 * Explore the tree by following the sequence.
 * Do not add anything the tree ("addSequenceAdd" does it).
 * Start at node cellId and character sequenceId.
 */
uint64_t addSequenceFollow (tree_t *tree, uint64_t cellId, char *binarySequence, int sequenceId) {
  assert(sequenceId >= 0);
  assert(cellId < tree->nCells);
  size_t edgeLength = 0;
  unsigned short sequenceNucleotide, edgeId, edgeNucleotide;
  edge_t *edge = NULL;
  //printf("Following sequence %s @%" PRIu64 "\n", sequence, cellId); fflush(stdout);
  for (; sequenceId >= 0; --sequenceId) {
    sequenceNucleotide = binarySequence[sequenceId];
    if (edge == NULL) {
      //printf("  new edge\n"); fflush(stdout);
      edgeId = sequenceNucleotide;
      edge = &tree->cells[cellId].edges[edgeId];
    }
    //printf("  seq id: %d, edge: %p (%i/%zu) -> %" PRIu64 ", nucleotide: %c\n", sequenceId, edge, edge->length, edgeLength, edge->cellId, DNA5_TO_CHAR[sequenceNucleotide]); fflush(stdout);
    assert(sequenceId >= 0);
    if (edge->cellId == NO_DATA) {
      return addSequenceAdd(tree, cellId, binarySequence, sequenceId, edge);
    }
    edgeNucleotide = getEdgeNucleotide(edge, edgeLength);
    //printf("  edge %p, value: %i, len: %zu/%u, nucleotide: %c (%u)\n", edge, edge->sequence, edgeLength, edge->length, DNA5_TO_CHAR[edgeNucleotide], edgeNucleotide); fflush(stdout);
    if (sequenceNucleotide == edgeNucleotide) {
      //printf("  following to %" PRIu64 "\n", edge->cellId); fflush(stdout);
      ++edgeLength;
      if (edgeLength == edge->length) {
        cellId = edge->cellId;
        //printf("  to the end -> %" PRIu64 "\n", cellId); fflush(stdout);
        edge = NULL;
        edgeLength = 0;
      }
    }
    else {
      //printf("  split @ %zu\n", edgeLength); fflush(stdout);
      if (edgeLength == 0) {
        cellId = edge->cellId;
        edge   = &tree->cells[cellId].edges[sequenceNucleotide];
        //printf("  edge is %p\n", edge); fflush(stdout);
        return addSequenceAdd(tree, cellId, binarySequence, sequenceId, edge);
      }
      else {
        cellId = splitEdgeTree(tree, cellId, edgeId, edgeLength, edgeNucleotide);
        edge   = &tree->cells[cellId].edges[sequenceNucleotide];
        ++tree->nEdges;
        return addSequenceAdd(tree, cellId, binarySequence, sequenceId, edge);
      }
    }
  }
  // set node
  //printf("  over with %zu\n", edgeLength); fflush(stdout);
  if (edgeLength != 0) {
    ++tree->nEdges;
    cellId = splitEdgeTree(tree, cellId, edgeId, edgeLength, N_NUCLEOTIDES);
  }
  return cellId;
}

bool addSequence (tree_t *tree, size_t l, char *sequence, char *quality, unsigned int fileId) {
  uint64_t cellId = 0;
  int sequenceId = l - 1;
  assert(strlen(sequence) == strlen(quality));
  assert(strlen(quality) == l);
  // Skip all the reads with a size lower than tree base.
  // Add one, because the tree base cells are numbered...
  if (l < TREE_BASE_SIZE + 1) {
    return false;
  }
  // convert to binary beforehand, because ambiguous nt are converted here
  convertToBinary(l, sequence);
  for (int i = 0; i < TREE_BASE_SIZE; ++i, --sequenceId) {
    cellId <<= NUCLEOTIDES_BITS;
    cellId += sequence[sequenceId];
  }
  assert(cellId < N_TREE_BASE);
  //printf("First id: %lu, %s\n", cellId, sequence);
  assert(cellId < N_TREE_BASE);
  cellId = addSequenceFollow(tree, cellId, sequence, sequenceId);
  setQuality(tree, cellId, l, quality, fileId);
  tree->depth = MAX(tree->depth, l);
  return true;
}

/**
 * Print the last nucleotides (close to the leaves) of the tree.
 */
void __printTree (const tree_t *tree, FILE *outFile, uint64_t *readId, char *read, size_t readPos, uint64_t cellId) {
  uint64_t nextCellId;
  edge_t *edge;
  cell_t *cell = &tree->cells[cellId];
  char *quality;
  //printf("Read: %s at %"PRIu64 "\n", read+tree->depth-readPos, cellId);
  if ((quality = findQuality(&tree->qualities, cellId)) != NULL) {
    //printf("\tGot quality\n"); fflush(stdout);
    ++(*readId);
    fprintf(outFile, "@read%" PRIu64 "_", *readId);
    for (unsigned int fileId = 0; fileId < parameters->nReadsFiles; ++fileId) {
      fprintf(outFile, "_%" PRIu32, cell->counts[fileId]);
    }
    fprintf(outFile, "\n%s\n+\n%s\n", read+tree->depth-readPos, quality);
    //assert(strlen(read+tree->depth-readPos) == strlen(cell->quality));
  }
  for (unsigned short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    edge = &cell->edges[nucleotide];
    nextCellId = edge->cellId;
    if (nextCellId != NO_DATA) {
      for (size_t i = 0; i < edge->length; ++i) {
        assert(1 + readPos + i <= tree->depth);
        read[tree->depth-readPos-1-i] = DNA5_TO_CHAR[getEdgeNucleotide(edge, i)];
      }
      __printTree(tree, outFile, readId, read, readPos + edge->length, nextCellId);
    }
  }
}

/**
 * Print the first nucleotides (close to the root) of the tree, then call __printTree.
 */
void _printTree (const tree_t *tree, FILE *outFile, uint64_t *readId, char *read, size_t readPos, uint64_t cellId) {
  if (readPos == TREE_BASE_SIZE) {
    __printTree (tree, outFile, readId, read, readPos, cellId);
    return;
  }
  assert(cellId < N_TREE_BASE);
  cellId <<= NUCLEOTIDES_BITS;
  for (unsigned short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    read[tree->depth-readPos-1] = DNA5_TO_CHAR[nucleotide];
    _printTree(tree, outFile, readId, read, readPos+1, cellId+nucleotide);
  }
}

/**
 * Open/close file, allocate the memory, and call _printTree
 */
void printTree (char *fileName, const tree_t *tree) {
  FILE *outFile = fopen(fileName, "w");
  if (outFile == NULL) {
    exit(EXIT_FAILURE);
  }
  uint64_t readId = 0;
  char *read = (char *) malloc((tree->depth+1) * sizeof(char));
  read[tree->depth] = 0;
  _printTree(tree, outFile, &readId, read, 0, 0);
  free(read);
  fclose(outFile);
  printf("Done with print tree.\n");
}

/**
 * Count the 3-mer, and possibly filter
 */
bool __filterTree (const tree_t *tree, size_t readPos, uint64_t cellId, unsigned short prevTriplet, count_t *prevTripletCount) {
  edge_t *edge;
  uint64_t nextCellId;
  bool foundRead = false;
  bool thresholdReached = false;
  cell_t *cell = &tree->cells[cellId];
  count_t tripletCount [N_TRIPLETS];
  unsigned short triplet;
  size_t prevReadPos = readPos;
  unsigned short nucleotide;
  //printf("\tCurrent triplet: %u @ readPos: %zu\n", prevTriplet, readPos);
  if (findQuality(&tree->qualities, cellId) != NULL) {
    //printf("\t\tHas quality\n");
    foundRead = true;
  }
  for (unsigned short edgeId = 0; edgeId < N_NUCLEOTIDES; ++edgeId) {
    edge = &cell->edges[edgeId];
    nextCellId = edge->cellId;
    //printf("\t\tTrying nucleotide %i\n", edgeId);
    if (nextCellId != NO_DATA) {
      thresholdReached = false;
      triplet = prevTriplet;
      memcpy(tripletCount, prevTripletCount, N_TRIPLETS * sizeof(count_t));
      readPos = prevReadPos;
      for (size_t edgeLength = 0; (edgeLength < edge->length) && (! thresholdReached); ++edgeLength, ++readPos) {
        nucleotide = getEdgeNucleotide(edge, edgeLength);
        triplet = ((triplet & TRIPLET_MASK) << NUCLEOTIDES_BITS) | nucleotide;
        //printf("\t\t\ttriplet: %u @ edgelen: %zu, nucleotide: %i\n", triplet, edgeLength, nucleotide);
        if (readPos >= TRIPLET-1) tripletCount[triplet] += 1;
        if (tripletCount[triplet] > parameters->lowComplexityThreshold) {
          //printf("\t\t\t\tDirect filtering\n");
          thresholdReached = true;
        }
     }
     if (thresholdReached) {
        unsetEdge(edge);
      }
      else {
        if (! __filterTree(tree, readPos+1, nextCellId, triplet, tripletCount)) {
          //printf("\t\tSecond filtering @ %zu\n", readPos);
          unsetEdge(edge);
        }
        else {
          foundRead = true;
        }
      }
    }
  }
  return foundRead;
}

/**
 * Initialize counts and run _filterTree
 */
bool _filterTree (const tree_t *tree, size_t readPos, uint64_t cellId, unsigned short triplet, count_t *tripletCount) {
  if (readPos == TREE_BASE_SIZE) {
    return __filterTree(tree, readPos, cellId, triplet, tripletCount);
  }
  unsigned short nextTriplet;
  bool foundRead = false;
  cellId <<= NUCLEOTIDES_BITS;
  triplet &= TRIPLET_MASK;
  triplet <<= NUCLEOTIDES_BITS;
  for (unsigned short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    nextTriplet = triplet | nucleotide;
    //printf("\tFirst filter @ %zu with nt %i and triplet %i\n", readPos, nucleotide, triplet);
    if (readPos >= TRIPLET-1) tripletCount[nextTriplet] += 1;
    if (_filterTree(tree, readPos+1, cellId | nucleotide, nextTriplet, tripletCount)) {
      foundRead = true;
    }
    if (readPos >= TRIPLET-1) tripletCount[nextTriplet] -= 1;
  }
  return foundRead;
}

/**
 * Initialize counts and run _filterTree
 */
unsigned int filterTree (const tree_t *tree) {
  count_t tripletCount [N_TRIPLETS];
  for (unsigned int i = 0; i < N_TRIPLETS; ++i) tripletCount[i] = 0;
  return _filterTree(tree, 0, 0, 0, tripletCount);
}

/**
 * Alter the edges so that they are comparable: either one is empty, or they should be equal.
 */
void mergeTreeEdge (tree_t *tree1, uint64_t cellId1, unsigned short edgeId1, tree_t *tree2, uint64_t cellId2, unsigned short edgeId2) {
  //printf("Merging edge ");
  //printEdge(edge1);
  //printf(" with ");
  //printEdge(edge2);
  //printf("\n"); fflush(stdout);
  unsigned int sequenceId;
  edge_t *edge1 = &tree1->cells[cellId1].edges[edgeId1];
  edge_t *edge2 = &tree2->cells[cellId2].edges[edgeId2];
  // first case: the first edge is empty
  if (! isSetEdge(edge1)) {
    ++tree1->nEdges;
    edge1->length   = edge2->length;
    edge1->sequence = edge2->sequence;
    edge1->cellId   = tree1->nCells;
    // Adding cells may trigger reallocation and invalidate edge pointer
    addCell(tree1);
    //printf("First case: ");
    //printEdge(edge1);
    //printf("\n"); fflush(stdout);
    return;
  }
  // go the first nucleotide which differs
  for (sequenceId = 0; (sequenceId < MIN(edge1->length, edge2->length)) && (getEdgeNucleotide(edge1, sequenceId) == getEdgeNucleotide(edge2, sequenceId)); ++sequenceId) {
    // do nothing
  }
  assert(sequenceId != 0);
  //printf("break point is %u\n", sequenceId); fflush(stdout);
  // split the edges
  if (sequenceId != edge1->length) {
    ++tree1->nEdges;
    splitEdgeTree(tree1, cellId1, edgeId1, sequenceId, N_NUCLEOTIDES);
    //printf("splitting edge 1: ");
    //printEdge(edge1);
    //printf("\n");
    //printCell(&tree1->cells[edge1->cellId]);
    //printf("\n"); fflush(stdout);
  }
  if (sequenceId != edge2->length) {
    splitEdgeTree(tree2, cellId2, edgeId2, sequenceId, N_NUCLEOTIDES);
    //printf("splitting edge 2: "); fflush(stdout);
    //printEdge(edge2);
    //printf("\n");
    //printCell(&tree2->cells[edge2->cellId]);
    //printf("\n"); fflush(stdout);
  }
}

/**
 * Add the counts of the second tree to the first.
 * Merge the qualities.
 */
void mergeTreeNodes (tree_t *tree1, uint64_t cellId1, tree_t *tree2, uint64_t cellId2) {
  cell_t *cell1 = &tree1->cells[cellId1];
  cell_t *cell2 = &tree2->cells[cellId2];
  char *quality1 = findQuality(&tree1->qualities, cellId1);
  char *quality2 = findQuality(&tree2->qualities, cellId2);
  for (size_t countId = 0; countId < parameters->nReadsFiles; ++countId) {
    cell1->counts[countId] += cell2->counts[countId];
  }
  if (quality2 == NULL) {
    return;
  }
  if (quality1 == NULL) {
    addQuality(&tree1->qualities, cellId1, strlen(quality2), quality2);
    return;
  }
  mergeQualities(quality1, quality2, strlen(quality1));
}

/**
 * Add the second tree to the first.
 * Merge nodes and proceed to the children edges.
 */
void _mergeTree (tree_t *tree1, uint64_t cellId1, tree_t *tree2, uint64_t cellId2) {
  //cell_t *cell1 = &tree1->cells[cellId1];
  //cell_t *cell2 = &tree2->cells[cellId2];
  //edge_t *edge1, *edge2;
  //printf("Merge tree %"PRIu64" (%p) vs %"PRIu64" (%p)\n", cellId1, cell1, cellId2, cell2); fflush(stdout);
  //printf("Node 1: ");
  //printCell(cell1);
  //printf("\nNode 2: ");
  //printCell(cell2);
  //printf("\n"); fflush(stdout);
  mergeTreeNodes(tree1, cellId1, tree2, cellId2);
  for (unsigned short nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    //edge2 = &cell2->edges[nucleotide];
    if (isSetEdge(&tree2->cells[cellId2].edges[nucleotide])) {
      //edge1 = &cell1->edges[nucleotide];
      // Adding node can trigger reallocation and invalidate edge pointer
      mergeTreeEdge(tree1, cellId1, nucleotide, tree2, cellId2, nucleotide);
      _mergeTree(tree1, tree1->cells[cellId1].edges[nucleotide].cellId, tree2, tree2->cells[cellId2].edges[nucleotide].cellId);
    }
  }
}

/**
 * Add the second tree to the first
 * Here, simply do it on all the sub-trees.
 */
void mergeTree (tree_t *tree1, tree_t *tree2) {
  for (uint64_t i = 0; i < N_TREE_BASE; ++i) {
    _mergeTree(tree1, i, tree2, i);
  }
  tree1->depth = MAX(tree1->depth, tree2->depth);
}

bwaidx_t *loadGenomeFile (char *indexName) {
  bwaidx_t *idx = bwa_idx_load(indexName, BWA_IDX_ALL);
  if (idx == NULL) {
    fprintf(stderr, "Index load failed.\n");
  }
  return idx;
}

#endif
