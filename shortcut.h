#ifndef SHORTCUT_H
#define SHORTCUT_H

#include "constants.h"

typedef struct {
  bool         isSet;
  size_t       depth;
  int          isRev;
  int          rid;
  unsigned int nMisses;
  bwtint_t     pos;
} shortCut_t;

void unsetShortCut (shortCut_t *shortCut) {
  shortCut->depth   = DEPTH_SHORT_CUT;
  shortCut->isSet   = false;
  shortCut->nMisses = 0;
}

shortCut_t *initializeShortCut () {
  shortCut_t *shortCut = (shortCut_t *) malloc(sizeof(shortCut_t));
  unsetShortCut(shortCut);
  return shortCut;
}

void setShortCut (shortCut_t *shortCut, bwtint_t k, size_t depth) {
  shortCut->isSet = true;
  shortCut->depth = depth;
  shortCut->pos   = bwt_sa(bwt, k);
  shortCut->pos   = bns_depos(bns, shortCut->pos, &shortCut->isRev);
  shortCut->rid   = bns_pos2rid(bns, shortCut->pos);
}

void resetShortCut (shortCut_t *shortCut, size_t depth) {
  shortCut->pos = (shortCut->isRev)? shortCut->pos+(depth-shortCut->depth): shortCut->pos-(depth-shortCut->depth);
  shortCut->depth = depth;
}

void incShortCut (shortCut_t *shortCut) {
  shortCut->pos = (shortCut->isRev)? shortCut->pos+1: shortCut->pos-1;
  ++shortCut->depth;
}

void addMiss (shortCut_t *shortCut) {
  shortCut->depth += shortCut->nMisses;
  ++shortCut->nMisses;
}

#endif
