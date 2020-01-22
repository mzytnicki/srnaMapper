// ==========================================================================
//                              srnaMapper
// ==========================================================================
// Copyright (C) 2019 Matthias Zytnicki, INRA
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Matthias Zytnicki or the INRA nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL MATTHIAS ZYTNICKI NOTR THE INRA BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Matthias Zytnicki <matthias.zytnicki@inra.fr>
// ==========================================================================

#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include <locale.h>

#include "constants.h"
#include "helper.h"
#include "stats.h"
#include "parameters.h"
#include "sam.h"
#include "edge.h"
#include "cell.h"
#include "quality.h"
#include "tree.h"
#include "bwt.h"
#include "state.h"
#include "sw.h"
#include "states.h"
#include "shortcut.h"
#include "path.h"
#include "main.h"



int main(int argc, char const ** argv) {
  int returnCode = 0;
  parameters_t param;
  stats_t stat;
  tree_t tree;
  tree2_t tree2;
  outputSam_t outputSam;
  bwaidx_t *idx = NULL;
  nReads = 0;
  parameters = &param;
  stats = &stat;
  initializeStats();
  returnCode = parseCommandLine(argc, argv);
  if (returnCode != EXIT_SUCCESS) return returnCode;
  createTree(&tree);
  for (unsigned int fileId = 0; fileId < parameters->nReadsFiles; ++fileId) {
    printf("Reading reads file %s...\n", parameters->readsFileNames[fileId]);
    returnCode = readReadsFile(parameters->readsFileNames[fileId], &tree, fileId);
    if (returnCode != EXIT_SUCCESS) return returnCode;
    puts("... done.");
  }
  puts("Filtering tree...");
  filterTree(&tree);
  printf("... done.\n");
  printf("Maximum read size: %zu\n", tree.depth);
  computeTreeStats(&tree);
  if (parameters->outputReadsFileName != NULL) {
    puts("Printing tree...");
    returnCode = printTree(parameters->outputReadsFileName, &tree);
    if (returnCode != EXIT_SUCCESS) return returnCode;
    puts("... done.");
  }
  puts("Loading genome...");
  puts("Simplifying tree...");
  copyTree(&tree2, &tree);
  freeTree(&tree);
  puts("... done");
  idx = loadGenomeFile(parameters->genomeFileName);
  if (idx == NULL) return EXIT_FAILURE;
  puts("... done.");
  puts("Mapping...");
  pac = idx->pac;
  bwt = idx->bwt;
  bns = idx->bns;
  FILE *outputSamFile = openSamFile();
  if (outputSamFile == NULL) {
    printf("Error!  Cannot write to output SAM file '%s'.\nExiting.\n", param.outputSamFileName);
  }
  outputSam.file = outputSamFile;
  createOutputSam(&outputSam, tree2.depth);
  map(&tree2, &outputSam);
  freeTree2(&tree2);
  freeOutputSam(&outputSam);
  bwa_idx_destroy(idx);
  puts("... done.");
  printStats();
  return 0;
}
