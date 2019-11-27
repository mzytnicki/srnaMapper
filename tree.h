#ifndef TREE_H
#define TREE_H

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
  uint64_t  nAllocated;
  quality_t qualities;
} tree_t;

void createTree (tree_t *tree) {
  tree->nAllocated = INIT_N_CELLS;
  tree->cells = (cell_t *) malloc(tree->nAllocated * sizeof(cell_t));
  tree->depth = TREE_BASE_SIZE;
  tree->nCells = N_TREE_BASE;
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
  tree->nAllocated = 0;
  tree->depth = 0;
  tree->nCells = 0;
  for (uint64_t i = 0; i < tree->nCells; ++i) {
    freeCell(&tree->cells[i]);
  }
  free(tree->cells);
  freeQualities(&tree->qualities);
}

void setQuality (tree_t *tree, size_t cellId, size_t l, char *quality, unsigned int fileId) {
  ++tree->cells[cellId].counts[fileId];
  _setQuality(&tree->qualities, cellId, l, quality);
}

uint64_t addCell (tree_t *tree) {
  if (tree->nCells == tree->nAllocated) {
    tree->nAllocated *= 2;
    //printf("reallocating cells...\n");
    if ((tree->cells = (cell_t *) realloc(tree->cells, tree->nAllocated * sizeof(cell_t))) == NULL) {
      printf("Cannot allocate memory for tree of size %lu.\nExiting.\n", tree->nAllocated);
      exit(EXIT_FAILURE);
    }
  }
  //printf("ncells: %zu / %zu\n", tree->nCells, tree->nAllocated);
  createCell(&tree->cells[tree->nCells]);
  ++tree->nCells;
  return tree->nCells-1;
}

/**
 * Create a new node, split an edge into two, and add the remaining part of the edge to then new cell.
 */
uint64_t splitEdgeTree (tree_t *tree, edge_t *edge, sequence_t length, unsigned short newEdgeId) {
  uint64_t newCellId = addCell(tree);
  splitEdge(edge, &tree->cells[newCellId], newCellId, length, newEdgeId);
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
 * Add the rest of the sequence and extend the edge
 */
uint64_t addSequenceAdd (tree_t *tree, uint64_t cellId, char *sequence, int sequenceId, edge_t *edge) {
  unsigned short nucleotide;
  //printf("Adding sequence %s @%" PRIu64 "\n", sequence, cellId);
  for (; sequenceId >= 0; --sequenceId) {
    nucleotide = CHAR_TO_DNA5[(int) sequence[sequenceId]];
    if (edge == NULL) {
      edge = &tree->cells[cellId].edges[nucleotide];
    }
    //printf("  seq id: %d, edge: %p (%i), nucleotide: %c\n", sequenceId, edge, edge->length, DNA5_TO_CHAR[nucleotide]);
    addEdgeNucleotide(edge, nucleotide);
    if (edge->length == MAX_EDGE_LENGTH) {
      cellId = addCell(tree);
      //printf("  to the end: %" PRIu64 "\n", cellId);
      setCellIdEdge(edge, cellId);
      edge = NULL;
    }
  }
  // set node
  if (edge != NULL) {
      cellId = addCell(tree);
      setCellIdEdge(edge, cellId);
  }
  //printf("  over with %" PRIu64 "\n", cellId);
  return cellId;
  /*
  cell_t *cell = &tree->cells[cellId];
  uint64_t newCellId = cell->children[childId];
  if (newCellId != NO_DATA) {
    return newCellId;
  }
  newCellId = addCell(tree);
  cell = &tree->cells[cellId]; // may be reallocated!
  cell->children[childId] = newCellId;
  return newCellId;
  */
}


uint64_t addSequenceFollow (tree_t *tree, uint64_t cellId, char *sequence, int sequenceId) {
  size_t edgeLength = 0;
  unsigned short sequenceNucleotide, edgeNucleotide;
  edge_t *edge = NULL;
  //printf("Following sequence %s @%" PRIu64 "\n", sequence, cellId); fflush(stdout);
  for (; sequenceId >= 0; --sequenceId) {
    sequenceNucleotide = CHAR_TO_DNA5[(int) sequence[sequenceId]];
    if (edge == NULL) {
      //printf("  new edge\n"); fflush(stdout);
      edge = &tree->cells[cellId].edges[sequenceNucleotide];
    }
    //printf("  seq id: %d, edge: %p (%i/%zu) -> %" PRIu64 ", nucleotide: %c\n", sequenceId, edge, edge->length, edgeLength, edge->cellId, DNA5_TO_CHAR[sequenceNucleotide]); fflush(stdout);
    if (edge->cellId == NO_DATA) {
      return addSequenceAdd(tree, cellId, sequence, sequenceId, edge);
    }
    edgeNucleotide = getEdgeNucleotide(edge, edgeLength); fflush(stdout);
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
      cellId = edge->cellId;
      if (edgeLength == 0) {
        edge = &tree->cells[cellId].edges[sequenceNucleotide];
        //printf("  edge is %p\n", edge); fflush(stdout);
        return addSequenceAdd(tree, cellId, sequence, sequenceId, edge);
      }
      else {
        cellId = splitEdgeTree(tree, edge, edgeLength, edgeNucleotide);
        edge   = &tree->cells[cellId].edges[sequenceNucleotide];
        return addSequenceAdd(tree, cellId, sequence, sequenceId, edge);
      }
    }
  }
  // set node
  //printf("  over with %zu\n", edgeLength); fflush(stdout);
  if (edgeLength != 0) {
    cellId = splitEdgeTree(tree, edge, edgeLength, N_NUCLEOTIDES);
  }
  return cellId;
}

bool addSequence (tree_t *tree, size_t l, char *sequence, char *quality, unsigned int fileId) {
  uint64_t cellId = 0;
  int sequenceId = l - 1;
  assert(strlen(sequence) == strlen(quality));
  assert(strlen(quality) == l);
  if (l < TREE_BASE_SIZE) {
    return false;
  }
  for (int i = 0; i < TREE_BASE_SIZE; ++i, --sequenceId) {
    cellId <<= NUCLEOTIDES_BITS;
    cellId += CHAR_TO_DNA5[(int) sequence[sequenceId]];
  }
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
      fprintf(outFile, "_%lu", cell->counts[fileId]);
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
int printTree (char *fileName, const tree_t *tree) {
  FILE *outFile = fopen(fileName, "w");
  if (outFile == NULL) return EXIT_FAILURE;
  uint64_t readId = 0;
  char *read = (char *) malloc((tree->depth+1) * sizeof(char));
  read[tree->depth] = 0;
  _printTree(tree, outFile, &readId, read, 0, 0);
  free(read);
  fclose(outFile);
  printf("Done with print tree.\n");
  return EXIT_SUCCESS;
}

int readReadsFile (char *fileName, tree_t *tree, unsigned int fileId) {
  FILE *inFile;
  char *line = NULL;
  char *sequence = NULL;
  char *quality = NULL;
  size_t len = 0;
  ssize_t nRead;
  inFile = fopen(fileName, "r");
  if (inFile == NULL) return EXIT_FAILURE;
  while ((nRead = getline(&line, &len, inFile)) != -1) {
    nRead = getline(&sequence, &len, inFile);
    if (nRead == -1) {
      fprintf(stderr, "Input file '%s' is corrupted.\nAborting.\n", fileName);
      return EXIT_FAILURE;
    }
    nRead = getline(&line, &len, inFile);
    if (nRead == -1) {
      fprintf(stderr, "Input file '%s' is corrupted.\nAborting.\n", fileName);
      return EXIT_FAILURE;
    }
    nRead = getline(&quality, &len, inFile);
    if (nRead == -1) {
      fprintf(stderr, "Input file '%s' is corrupted.\nAborting.\n", fileName);
      return EXIT_FAILURE;
    }
    assert(strlen(sequence) == strlen(quality));
    assert(strlen(sequence) == (unsigned long) nRead);
    trimSequence(nRead, sequence);
    trimSequence(nRead, quality);
    ++stats->nReads;
    if (! addSequence(tree, nRead-1, sequence, quality, fileId)) {
      ++stats->nShortReads;
    }
  }
  free(line);
  free(sequence);
  free(quality);
  fclose(inFile);
  return EXIT_SUCCESS;
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

bwaidx_t *loadGenomeFile (char *indexName) {
  bwaidx_t *idx = bwa_idx_load(indexName, BWA_IDX_ALL);
  if (idx == NULL) {
    fprintf(stderr, "Index load failed.\n");
  }
  return idx;
}

#endif
