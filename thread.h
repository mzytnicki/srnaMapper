#ifndef THREAD_H
#define THREAD_H

#include "constants.h"
#include "parameters.h"
#include "tree2.h"
#include "main.h"

#define THREAD_CELLID_STEP 0x100

/**
 * Data which are given to the threads
 * tree:        the reads tree
 * samFile:     the name of the output sam file
 * firstCellId: the id of the first cell scanned by this thread at depth TREE_BASE_SIZE (included)
 * lastCellId:  the id of the last  cell scanned by this thread at depth TREE_BASE_SIZE (included)
 */
typedef struct {
  tree2_t         *tree;
  FILE            *samFile;
  pthread_mutex_t *incMutex;
  pthread_mutex_t *writeMutex;
} threadParameters_t;

typedef struct {
  pthread_t *threads;
} thread_t;

unsigned int threadStep;

void initializeThreads (thread_t *threads) {
  threads->threads = (pthread_t *) malloc((parameters->nThreads-1) * sizeof(pthread_t));
}

void freeThreads (thread_t *threads) {
  free(threads->threads);
}

void updateBounds (uint32_t *firstCellId, uint32_t *lastCellId, pthread_mutex_t *incMutex) {
  //printf("Acquiring i mutex\n"); fflush(stdout);
  pthread_mutex_lock(incMutex);
  //printf("Acquiring i mutex done\n"); fflush(stdout);
  *firstCellId = threadStep * THREAD_CELLID_STEP;
  *lastCellId  = (threadStep + 1) * THREAD_CELLID_STEP - 1;
  ++threadStep;
  //printf("Releasing i mutex\n"); fflush(stdout);
  pthread_mutex_unlock(incMutex);
  //printf("Releasing i mutex done\n"); fflush(stdout);
}

void *threadMain (void *parametersVoid) {
  //printf("Preparing thread...\n"); fflush(stdout);
  threadParameters_t *parameters  = (threadParameters_t *) parametersVoid;
  states_t           *states      = initializeStates(parameters->tree->depth);
  path_t             *path        = initializePath(parameters->tree->depth);
  uint32_t            firstCellId;
  uint32_t            lastCellId;
  cellVisitor_t cellVisitor;
  outputSam_t   outputSam;
  outputSam.file = parameters->samFile;
  createOutputSam(&outputSam, parameters->tree->depth, parameters->writeMutex);
  updateBounds(&firstCellId, &lastCellId, parameters->incMutex);
  clearCellVisitor(&cellVisitor, &parameters->tree->cellInfos);
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
    writeToSam(&outputSam);
    //printf("Ending thread %" PRIu32 "-%" PRIu32 " at depth %zu with cell #%" PRIu32 "\n", parameters->firstCellId, parameters->lastCellId, path->nCells, path->cellIds[path->nCells]); fflush(stdout);
    updateBounds(&firstCellId, &lastCellId, parameters->incMutex);
  }
  freeOutputSam(&outputSam);
  freeStates(states);
  freePath(path);
  return NULL;
}

void startThreads (tree2_t *tree, FILE *samFile) {
  pthread_mutex_t incMutex   = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_t writeMutex = PTHREAD_MUTEX_INITIALIZER;
  thread_t threads;
  threadParameters_t threadParameters;
  if (pthread_mutex_init(&incMutex, NULL) != 0) {
    fprintf(stderr, "Error! Cannot initialize mutex.\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  if (pthread_mutex_init(&writeMutex, NULL) != 0) {
    fprintf(stderr, "Error! Cannot initialize mutex.\nExiting.\n");
    exit(EXIT_FAILURE);
  }
  initializeThreads(&threads);
  threadParameters.tree       =  tree;
  threadParameters.samFile    =  samFile;
  threadParameters.incMutex   = &incMutex;
  threadParameters.writeMutex = &writeMutex;
  threadStep = 0;
  //printf("Got %u threads\n", parameters->nThreads); fflush(stdout);
  for (unsigned int threadId = 0; threadId < parameters->nThreads-1; ++threadId) {
    if (pthread_create(&threads.threads[threadId], NULL, threadMain, (void *) &threadParameters) != 0) {
      fprintf(stderr, "Error!  Thread %u cannot be created.\nExiting.\n", threadId);
      exit(EXIT_FAILURE);
    }
  }
  // Use the main thread too
  threadMain((void *) &threadParameters);
  for (unsigned int threadId = 0; threadId < parameters->nThreads-1; ++threadId) {
    pthread_join(threads.threads[threadId], NULL); 
  }
  //pthread_exit(NULL);
  pthread_mutex_destroy(&incMutex);
  pthread_mutex_destroy(&writeMutex);
  freeThreads(&threads);
}

#endif
