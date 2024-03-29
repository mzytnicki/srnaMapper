#ifndef PATH_T
#define PATH_T

#include "constants.h"
#include "edge.h"
#include "cell.h"
#include "tree.h"
#include "tree2.h"
#include "states.h"

/******* Path type *******/
/**
 * The path stores the path to a cell of the reads tree.
 * It is:
 *   - the nucleotide chose up to that cell
 *   - the corresponding cellId for each cell of the path
 *   - the tree depth
 *   - the corresponding read
 *   - the index of the next char to be written for the read (not the one you want to change when you change the branch)
 *   - nCells is actually the index of the last used cellId
 */

typedef struct {
  unsigned short *nucleotides;
  unsigned short *children;
  uint32_t       *cellIds;
  size_t          nCells;
  size_t          maxDepth;
  size_t          depth;
  edge_t         *edges;
  size_t          edgeLength;
  char           *read;
  size_t          readPos;
} path_t;

void printPath (path_t *path) {
  assert(path->depth <= path->maxDepth);
  printf("Read: %s\n", path->read + path->readPos);
  for (size_t i = 0; i < path->depth; ++i) {
    assert(path->nucleotides[i] < N_NUCLEOTIDES);
    printf("%i ", path->nucleotides[i]);
  }
  printf("\n");
  for (size_t i = 0; i <= path->nCells; ++i) {
    printf("%" PRIu32 " ", path->cellIds[i]);
    //if ((i >= TREE_BASE_SIZE) && (i < path->nCells)) {
    if (i >= TREE_BASE_SIZE) {
      printf("child #%i ", path->children[i]);
      //printEdge2(&path->edges[i]);
      printEdge(&path->edges[i]);
    }
    printf("  ");
  }
  printf("\n(depth: %zu, # cells: %zu, read pos: %zu, cellId: %" PRIu32 ", edge length: %zu)\n", path->depth, path->nCells, path->readPos, path->cellIds[path->nCells], path->edgeLength); fflush(stdout);
}

void clearPath (path_t *path) {
  path->nCells         = 0;
  path->depth          = 0;
  path->edgeLength     = 0;
  path->readPos        = path->maxDepth;
}

path_t *initializePath (size_t maxDepth) {
  path_t *path         = (path_t *)         mallocOrDie(sizeof(path_t));
  path->nucleotides    = (unsigned short *) mallocOrDie(maxDepth * sizeof(unsigned short));
  path->children       = (unsigned short *) mallocOrDie(maxDepth * sizeof(unsigned short));
  path->cellIds        = (uint32_t *)       mallocOrDie((maxDepth+1) * sizeof(uint32_t));
  path->edges          = (edge_t *)         mallocOrDie((maxDepth+1) * sizeof(edge_t));
  //path->edges          = (edge2_t *)        malloc((maxDepth+1) * sizeof(edge2_t));
  path->read           = (char *)           mallocOrDie((maxDepth+1) * sizeof(char));
  path->maxDepth       = maxDepth;
  path->cellIds[0]     = 0;
  path->read[maxDepth] = 0;
  clearPath(path);
  return path;
}

void freePath (path_t *path) {
  free(path->nucleotides);
  free(path->children);
  free(path->cellIds);
  free(path->edges);
  free(path->read);
  free(path);
}

void appendNucleotidePath (path_t *path, unsigned short nucleotide, char c) {
  assert(path->depth < path->maxDepth);
  path->nucleotides[path->depth] = nucleotide;
  ++path->depth;
  --path->readPos;
  path->read[path->readPos] = c;
}

bool isAfterFirstCellId (size_t depth, size_t cellId, uint32_t firstCellId) {
  /*
  if ((firstCellId >> (NUCLEOTIDES_BITS * (TREE_BASE_SIZE - path->nCells))) > path->cellIds[path->nCells]) {
    printf("Depth: %zu, cell: %" PRIu32 ", 1st cell: %" PRIu32 ", computed: %" PRIu32 ", target sequence: ", path->nCells, path->cellIds[path->depth], firstCellId, firstCellId >> (NUCLEOTIDES_BITS * (TREE_BASE_SIZE - path->nCells)));
    printSequence(firstCellId, TREE_BASE_SIZE);
    printf(", current sequence: ");
    printSequence(firstCellId, path->nCells);
    printf(", computed sequence: ");
    printSequence(firstCellId >> (NUCLEOTIDES_BITS * (TREE_BASE_SIZE - path->nCells)), path->nCells);
    printf("\n");
    fflush(stdout);
  }
  */
  if (depth == 0) {
    return true;
  }
  if (depth > TREE_BASE_SIZE) {
    return true;
  }
  return ((firstCellId >> (NUCLEOTIDES_BITS * (TREE_BASE_SIZE - depth))) <= cellId);
}

/**
 * Step into the last cells of the tree in a DFS fashion.
 */
bool goDownTreeBase (path_t *path, uint32_t firstCellId) {
  uint32_t cellId = path->cellIds[path->nCells] << NUCLEOTIDES_BITS;
  for (unsigned char nucleotide = 0; nucleotide < N_NUCLEOTIDES; ++nucleotide) {
    if (isAfterFirstCellId(path->depth+1, cellId, firstCellId)) {
      path->cellIds[path->nCells+1] = cellId;
      appendNucleotidePath(path, nucleotide, "ACGT"[nucleotide]);
      ++path->nCells;
      return true;
    }
    ++cellId;
  }
  //printf("\t\tGoing to base\n");
  //printPath(path);
  return false;
}

/**
 * Step into the first cells of the tree in a DFS fashion.
 */
bool goDownTree2NotBase (const tree2_t *tree, path_t *path) {
  //printf("\t\tGoing to tree\n");
  //edge2_t *edge = NULL;
  edge_t *edge = NULL;
  unsigned short nucleotide;
  //printPath(path);
  assert(path->depth <= tree->depth);
  assert(path->nCells <= tree->depth);
  if (path->edgeLength == 0) {
    cell2_t *cell = &tree->cells[path->cellIds[path->nCells]];
    if (cell->nEdges == 0) {
      //printf("\t\t\tNothing\n");
      return false;
    }
    edge = &tree->edges[cell->firstEdge];
    //nucleotide = getEdge2Nucleotide(tree, edge, 0);
    nucleotide = getNucleotide(edge->sequence, 0);
    path->edges[path->nCells]    = *edge;
    path->children[path->nCells] = 0;
    path->edgeLength             = 1;
    appendNucleotidePath(path, nucleotide, DNA5_TO_CHAR[nucleotide]);
    if (path->edgeLength == edge->length) {
      path->edgeLength = 0;
      path->cellIds[++path->nCells] = edge->cellId;
    }
    //printf("to %zu with read '%s'\n", path->depth, path->read+path->readPos);
    return true;
  }
  edge = &path->edges[path->nCells];
  //nucleotide = getEdge2Nucleotide(tree, edge, path->edgeLength);
  nucleotide = getNucleotide(edge->sequence, path->edgeLength);
  ++path->edgeLength;
  appendNucleotidePath(path, nucleotide, DNA5_TO_CHAR[nucleotide]);
  if (path->edgeLength == edge->length) {
    path->edgeLength = 0;
    path->cellIds[++path->nCells] = edge->cellId;
  }
  return true;
}

/**
 * Step into the tree in a DFS fashion.
 * Return false if search is exhausted.
 */
bool goDownTree2 (const tree2_t *tree, path_t *path, uint32_t firstCellId) {
  //printf("\tGo down from height %zu with read '%s'\n", path->depth, path->read+path->readPos);
  //printPath(path);
  if (path->depth < TREE_BASE_SIZE) {
    return goDownTreeBase(path, firstCellId);
  }
  return goDownTree2NotBase(tree, path);
}

bool isBeforeLastCellId (path_t *path, uint32_t lastCellId) {
  /*
  if (path->cellIds[path->nCells] > (lastCellId >> (NUCLEOTIDES_BITS * (TREE_BASE_SIZE - path->nCells)))) {
    printf("Depth: %zu, cell: %" PRIu32 ", last cell: %" PRIu32 ", computed: %" PRIu32 ", sequence: ", path->nCells, path->cellIds[path->nCells], lastCellId, lastCellId >> (NUCLEOTIDES_BITS * (TREE_BASE_SIZE - path->nCells)));
    printSequence(lastCellId >> (NUCLEOTIDES_BITS * (TREE_BASE_SIZE - path->nCells)), path->nCells);
    printf("\n");
    fflush(stdout);
  }
  */
  if (path->depth == 0) {
    return true;
  }
  if (path->depth > TREE_BASE_SIZE) {
    return true;
  }
  return (path->cellIds[path->nCells] <= (lastCellId >> (NUCLEOTIDES_BITS * (TREE_BASE_SIZE - path->nCells))));
}


/**
 * Step into the last cells of the tree in a BFS fashion.
 */
bool goRightTreeBase (path_t *path, uint32_t lastCellId) {
  //uint64_t cellId = path->cellIds[path->nCells], nextCellId = cellId + 1, mask = NUCLEOTIDE_MASK;
  unsigned short newNucleotide = 0;
  assert(path->depth <= TREE_BASE_SIZE);
  assert(path->nCells <= TREE_BASE_SIZE);
  //printf("    Entering go right base will cell %" PRIu32 " at depth %zu, last nt %i, and read pos %zu\n", path->cellIds[path->nCells], path->depth, path->nucleotides[path->depth-1], path->readPos);
  for (; path->nucleotides[path->depth-1] == N_NUCLEOTIDES - 1; --path->depth, ++path->readPos) {
    if (path->depth == 1) return false;
  }
  newNucleotide = path->nucleotides[path->depth-1] + 1;
  path->read[path->readPos] = DNA5_TO_CHAR[newNucleotide];
  path->nucleotides[path->depth-1] = newNucleotide;
  /*
  for (unsigned int offset = 0; (cellId & mask) != (nextCellId & mask); ++offset, mask <<= NUCLEOTIDES_BITS) {
    newNucleotide = ((nextCellId & mask) >> (NUCLEOTIDES_BITS * offset));
    //printf("      Offset: %u, mask: %" PRIu64 ", depth: %zu, new char: %c\n", offset, mask, path->depth, DNA5_TO_CHAR[newNucleotide]);
    path->read[path->readPos+offset] = DNA5_TO_CHAR[newNucleotide];
    path->nucleotides[path->depth-offset-1] = newNucleotide;
  }
  */
  path->nCells = path->depth;
  path->cellIds[path->nCells] = path->cellIds[path->nCells] + 1;
  //printf("      Leaving base will cell %" PRIu32 ", nucleotide %c, read %s\n", path->cellIds[path->nCells], "ACGT"[newNucleotide], path->read + path->readPos);
  path->edgeLength = 0;
  return isBeforeLastCellId(path, lastCellId);
}

/**
 * Step into the first cells of the tree in a BFS fashion.
 */
bool goRightTree2NotBase (const tree2_t *tree, path_t *path, uint64_t lastCellId) {
  assert(path->depth <= path->maxDepth);
  assert(path->nCells <= path->maxDepth);
  //edge2_t *edge;
  edge_t *edge;
  unsigned short childId, nucleotide;
  //printf("    not base\n");
  while (true) {
    if (path->depth <= TREE_BASE_SIZE) {
      ++path->nCells;
      return goRightTreeBase(path, lastCellId);
    }
    //printf("    ... trying depth %zu, # cell %zu, edge len %zu\n", path->depth, path->nCells, path->edgeLength);
    assert(path->depth >= path->edgeLength);
    if (path->edgeLength == 0) {
      --path->nCells;
      path->edgeLength = path->edges[path->nCells].length;
    }
    path->depth   -= path->edgeLength;
    path->readPos += path->edgeLength;
    //printf("    ... trying depth %zu, # cell %zu, edge len %zu\n", path->depth, path->nCells, path->edgeLength);
    //printf("    ... cellId: %" PRIu32 "/%" PRIu32 ", nt: %i\n", path->cellIds[path->nCells], tree->nCells, path->nucleotides[path->depth]);
    //printf("    ");
    //printEdge(&tree->cells[path->cellIds[path->nCells]].edges[path->nucleotides[path->depth]]);
    assert(path->depth < path->maxDepth);
    // following assert does not work if mutation at first nucleotide...
    //assert(tree->cells[path->cellIds[path->nCells]].edges[path->nucleotides[path->depth]].cellId == path->edges[path->nCells].cellId);
    cell2_t *cell = &tree->cells[path->cellIds[path->nCells]];
    //printf("    ... cell: ");
    //printCell2(cell);
    //printf("\n"); fflush(stdout);
    childId = path->children[path->nCells] + 1;
    //printf("    ... child: %i/%i\n", childId, cell->nEdges);
    if (childId < cell->nEdges) {
      //printf("    ... trying nucleotide %c\n", "ACGTN"[nucleotide]);
      edge = &tree->edges[cell->firstEdge + childId];
      //printf("    ... edge ");
      //printEdge(edge);
      //printf("\n");
      //nucleotide = getEdge2Nucleotide(tree, edge, 0);
      nucleotide = getNucleotide(edge->sequence, 0);
      path->nucleotides[path->depth] = nucleotide;
      ++path->depth;
      --path->readPos;
      path->read[path->readPos] = DNA5_TO_CHAR[nucleotide];
      path->children[path->nCells] = childId;
      path->edges[path->nCells] = *edge;
      path->edgeLength = 1;
      if (path->edgeLength == edge->length) {
        path->edgeLength = 0;
        path->cellIds[++path->nCells] = edge->cellId;
      }
      //printf("\t\t\tFinally, path is\n");
      //printPath(path);
      return true;
    }
    --path->nCells;
    path->edgeLength = path->edges[path->nCells].length;
  }
  //printf("    ... going to nothing\n");
  return false;
}

/**
 * Step into the tree in a BFS fashion.
 * Return false if search is exhausted.
 */
bool goRightTree2 (const tree2_t *tree, path_t *path, uint32_t lastCellId) {
  //printf("  Go right read from depth %zu with read '%s'\n", path->depth, path->read+path->readPos);
  //printPath(path);
  if (path->depth > TREE_BASE_SIZE) {
    return goRightTree2NotBase(tree, path, lastCellId);
  }
  return goRightTreeBase(path, lastCellId);
}

/**
 * Go to next child node in reads tree.  If none, go to sibling
 */
bool goNextTree2 (const tree2_t *tree, states_t *states, path_t *path, uint32_t firstCellId, uint32_t lastCellId, bool mappable) {
  // First make sure that we are not before the first cell
  //printf("cell id %" PRIu32 "\n", path->cellIds[path->nCells]); fflush(stdout);
  //printf("Entering goNextTree2 with states %p, path %p, depth %zu, ncells %zu, cell id %" PRIu32 ", interval [%" PRIu32 "-%" PRIu32 "]\n", states, path, path->depth, path->nCells, path->cellIds[path->nCells], firstCellId, lastCellId); fflush(stdout);
  //printf("Entering goNext\n");
  //printPath(path);
  /*
  while (! isAfterFirstCellId(path, firstCellId)) {
    assert(path->nCells > 0);
    printf("  Goto next adjusted\n");
    if (! goRightTree2(tree, path, lastCellId)) {
      return false;
    }
    printPath(path);
  }
  */
  if (mappable) {
    //printf("    Mappable\n");
    if (goDownTree2(tree, path, firstCellId)) {
      //printf("    Going down\n");
      //printPath(path);
      return true;
    }
    if (path->depth == 0) {
      return false;
    }
  }
  //printf("  Going right\n");
  //printStates(states, path->depth);
  //printPath(path);
  if (! goRightTree2(tree, path, lastCellId)) {
    return false;
  }
  //printf("  Now\n");
  //printStates(states, path->depth);
  backtrackStates(states, path->depth);
  return true;
}

#endif
