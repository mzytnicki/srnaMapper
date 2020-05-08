#ifndef THREAD_H
#define THREAD_H

#include "constants.h"
#include "parameters.h"
#include "tree2.h"
#include "main.h"

#define THREAD_BYTES_STEP 0x10000
#define THREAD_CELLID_STEP 0x100

/**
 * Data which are given to the mapping threads
 * tree:        the reads tree
 * samFile:     the name of the output sam file
 * firstCellId: the id of the first cell scanned by this thread at depth TREE_BASE_SIZE (included)
 * lastCellId:  the id of the last  cell scanned by this thread at depth TREE_BASE_SIZE (included)
 */
typedef struct {
  tree_t          *tree;
  unsigned int     fileId;
  char            *fastqFileName;
  FILE            *fastqFile;
  pthread_mutex_t *incMutex;
} readingThreadParameters_t;

typedef struct {
  tree2_t         *tree;
  FILE           **samFiles;
  pthread_mutex_t *incMutex;
  pthread_mutex_t *writeMutex;
} mappingThreadParameters_t;

typedef struct {
  tree_t *tree1;
  tree_t *tree2;
} mergeThreadParameters_t;

typedef struct {
  pthread_t *threads;
} thread_t;

unsigned int threadStep;

void createThreads (thread_t *threads) {
  threads->threads = (pthread_t *) malloc((parameters->nThreads-1) * sizeof(pthread_t));
}

void freeThreads (thread_t *threads) {
  free(threads->threads);
}

void updateBounds (uint64_t *firstCellId, uint64_t *lastCellId, unsigned long step, pthread_mutex_t *incMutex) {
  pthread_mutex_lock(incMutex);
  *firstCellId = threadStep * step;
  *lastCellId  = (threadStep + 1) * step - 1;
  ++threadStep;
  pthread_mutex_unlock(incMutex);
}

void *mergingThreadMain (void *parametersVoid) {
  mergeThreadParameters_t *parameters  = (mergeThreadParameters_t *) parametersVoid;
  mergeTree(parameters->tree1, parameters->tree2);
  return NULL;
}

void startMergingThreads (thread_t *threads, tree_t *mainTree, tree_t *otherTrees) {
  unsigned int nTrees = parameters->nThreads;
  unsigned int nMerges = nTrees / 2;
  mergeThreadParameters_t *threadParameters = (mergeThreadParameters_t *) malloc(nMerges * sizeof(mergeThreadParameters_t));
  tree_t **allTrees = (tree_t **) malloc((parameters->nThreads) * sizeof(tree_t *));
  for (unsigned int threadId = 0; threadId < parameters->nThreads-1; ++threadId) {
    allTrees[threadId+1] = &otherTrees[threadId];
  }
  allTrees[0] = mainTree;
  while (nTrees != 1) {
    for (unsigned int threadId = 0; threadId < nMerges; ++threadId) {
      threadParameters[threadId].tree1 = allTrees[threadId*2];
      threadParameters[threadId].tree2 = allTrees[threadId*2+1];
      if (pthread_create(&threads->threads[threadId], NULL, mergingThreadMain, (void *) &threadParameters[threadId]) != 0) {
        fprintf(stderr, "Error!  Thread %u cannot be created.\nExiting.\n", threadId);
        exit(EXIT_FAILURE);
      }
    }
    for (unsigned int threadId = 0; threadId < nMerges; ++threadId) {
      pthread_join(threads->threads[threadId], NULL); 
    }
    nTrees = nMerges + (nTrees & 1);
    for (unsigned int threadId = 0; threadId < nTrees; ++threadId) {
      allTrees[threadId] = allTrees[2*threadId];
    }
    nMerges = nTrees / 2;
  }
  free(allTrees);
  free(threadParameters);
}

void *readingThreadMain (void *parametersVoid) {
  //printf("Preparing thread...\n"); fflush(stdout);
  readingThreadParameters_t *parameters  = (readingThreadParameters_t *) parametersVoid;
  uint64_t posStart, posEnd;
  bool stillData = true;
  updateBounds(&posStart, &posEnd, THREAD_BYTES_STEP, parameters->incMutex);
  while (stillData) {
    //printf("Starting thread %" PRIu32 "-%" PRIu32 ": ", firstCellId, lastCellId);
    //printSequence(firstCellId, TREE_BASE_SIZE);
    //printf("-");
    //printSequence(lastCellId, TREE_BASE_SIZE);
    //printf("\n");
    //fflush(stdout);
    stillData = readReadsFile(parameters->fastqFile, parameters->fastqFileName, parameters->tree, parameters->fileId, posStart, posEnd);
    //printf("Ending thread %" PRIu32 "-%" PRIu32 " at depth %zu with cell #%" PRIu32 "\n", parameters->firstCellId, parameters->lastCellId, path->nCells, path->cellIds[path->nCells]); fflush(stdout);
    updateBounds(&posStart, &posEnd, THREAD_BYTES_STEP, parameters->incMutex);
  }
  return NULL;
}

void startReadingThreads (thread_t *threads, tree_t *tree, FILE **fastqFiles) {
  pthread_mutex_t incMutex                    = PTHREAD_MUTEX_INITIALIZER;
  tree_t *trees                               = (tree_t *)                    malloc((parameters->nThreads-1) * sizeof(tree_t));
  readingThreadParameters_t *threadParameters = (readingThreadParameters_t *) malloc(parameters->nThreads     * sizeof(readingThreadParameters_t));
  if (pthread_mutex_init(&incMutex, NULL) != 0) {
    fprintf(stderr, "Error! Cannot initialize mutex.\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  for (unsigned int threadId = 0; threadId < parameters->nThreads-1; ++threadId) {
    createTree(&trees[threadId]);
    threadParameters[threadId].tree     = &trees[threadId];
    threadParameters[threadId].incMutex = &incMutex;
  }
  // last thread parameter is for the main thread
  threadParameters[parameters->nThreads-1].tree     = tree;
  threadParameters[parameters->nThreads-1].incMutex = &incMutex;
  //printf("Got %u threads\n", parameters->nThreads); fflush(stdout);
  for (size_t fileId = 0; fileId < parameters->nReadsFiles; ++fileId) {
    fprintf(stderr, "Reading file '%s'...\n", parameters->readsFileNames[fileId]);
    threadStep = 0;
    for (unsigned int threadId = 0; threadId < parameters->nThreads-1; ++threadId) {
      threadParameters[threadId].fastqFileName = parameters->readsFileNames[fileId];
      threadParameters[threadId].fastqFile     = fopen(parameters->readsFileNames[fileId], "r");
      threadParameters[threadId].fileId        = fileId;
      if (pthread_create(&threads->threads[threadId], NULL, readingThreadMain, (void *) &threadParameters[threadId]) != 0) {
        fprintf(stderr, "Error!  Thread %u cannot be created.\nExiting.\n", threadId);
        exit(EXIT_FAILURE);
      }
    }
    // Use the main thread too
    threadParameters[parameters->nThreads-1].fastqFileName = parameters->readsFileNames[fileId];
    threadParameters[parameters->nThreads-1].fastqFile     = fastqFiles[fileId];
    threadParameters[parameters->nThreads-1].fileId        = fileId;
    readingThreadMain((void *) &threadParameters[parameters->nThreads-1]);
    for (unsigned int threadId = 0; threadId < parameters->nThreads-1; ++threadId) {
      pthread_join(threads->threads[threadId], NULL); 
      fclose(threadParameters[threadId].fastqFile);
    }
    fprintf(stderr, "... done.\n");
  }
  // merge all trees to the main tree
  if (parameters->nThreads > 1) {
    startMergingThreads (threads, tree, trees);
  }
  for (unsigned int threadId = 0; threadId < parameters->nThreads-1; ++threadId) {
    freeTree(&trees[threadId]);
  }
  free(trees);
  free(threadParameters);
  pthread_mutex_destroy(&incMutex);
}

void *mappingThreadMain (void *parametersVoid) {
  //printf("Preparing thread...\n"); fflush(stdout);
  mappingThreadParameters_t *parameters  = (mappingThreadParameters_t *) parametersVoid;
  states_t                  *states      = initializeStates(parameters->tree->depth);
  path_t                    *path        = initializePath(parameters->tree->depth);
  uint64_t                   firstCellId;
  uint64_t                   lastCellId;
  cellVisitor_t              cellVisitor;
  outputSam_t                outputSam;
  outputSam.outputFiles = parameters->samFiles;
  createOutputSam(&outputSam, parameters->tree->depth, parameters->writeMutex);
  updateBounds(&firstCellId, &lastCellId, THREAD_CELLID_STEP, parameters->incMutex);
  clearCellVisitor(&cellVisitor);
  while (firstCellId <= N_TREE_BASE) {
    //printf("Starting thread %" PRIu32 "-%" PRIu32 ": ", firstCellId, lastCellId);
    //printSequence(firstCellId, TREE_BASE_SIZE);
    //printf("-");
    //printSequence(lastCellId, TREE_BASE_SIZE);
    //printf("\n");
    //fflush(stdout);
    //TODO: Should work without next line
    //clearCellVisitor(&cellVisitor, &parameters->tree->cellInfos);
    clearStates(states);
    clearPath(path);
    _map(parameters->tree, states, path, &cellVisitor, firstCellId, lastCellId, &outputSam);
    writeToSam(&outputSam, true);
    //printf("Ending thread %" PRIu32 "-%" PRIu32 " at depth %zu with cell #%" PRIu32 "\n", parameters->firstCellId, parameters->lastCellId, path->nCells, path->cellIds[path->nCells]); fflush(stdout);
    updateBounds(&firstCellId, &lastCellId, THREAD_CELLID_STEP, parameters->incMutex);
  }
  freeOutputSam(&outputSam);
  freeStates(states);
  freePath(path);
  return NULL;
}

void startMappingThreads (thread_t *threads, tree2_t *tree, FILE **samFiles) {
  pthread_mutex_t incMutex   = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_t writeMutex = PTHREAD_MUTEX_INITIALIZER;
  mappingThreadParameters_t threadParameters;
  if (pthread_mutex_init(&incMutex, NULL) != 0) {
    fprintf(stderr, "Error! Cannot initialize mutex.\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  if (pthread_mutex_init(&writeMutex, NULL) != 0) {
    fprintf(stderr, "Error! Cannot initialize mutex.\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  threadParameters.tree       =  tree;
  threadParameters.samFiles   =  samFiles;
  threadParameters.incMutex   = &incMutex;
  threadParameters.writeMutex = &writeMutex;
  threadStep = 0;
  //printf("Got %u threads\n", parameters->nThreads); fflush(stdout);
  for (unsigned int threadId = 0; threadId < parameters->nThreads-1; ++threadId) {
    if (pthread_create(&threads->threads[threadId], NULL, mappingThreadMain, (void *) &threadParameters) != 0) {
      fprintf(stderr, "Error!  Thread %u cannot be created.\nExiting.\n", threadId);
      exit(EXIT_FAILURE);
    }
  }
  // Use the main thread too
  mappingThreadMain((void *) &threadParameters);
  for (unsigned int threadId = 0; threadId < parameters->nThreads-1; ++threadId) {
    pthread_join(threads->threads[threadId], NULL); 
  }
  //pthread_exit(NULL);
  pthread_mutex_destroy(&incMutex);
  pthread_mutex_destroy(&writeMutex);
}

#endif
