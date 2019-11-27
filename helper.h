#ifndef HELPER_H
#define HELPER_H

#include "constants.h"

void trimSequence (size_t l, char *s) {
  s[l-1] = 0;
}

char *reverseSequence (char *sequence, char *revSequence, size_t size) {
  assert(strlen(sequence) == size);
  for (size_t i = 0; i < size; ++i) {
    revSequence[i] = sequence[size-i-1];
  }
  revSequence[size] = 0;
  return revSequence;
}

char *reverseComplementSequence (char *sequence, char *revSequence, size_t size) {
  for (size_t i = 0; i < size; ++i) {
    revSequence[i] = COMP[(int) sequence[size-i-1]];
  }
  revSequence[size] = 0;
  return revSequence;
}

unsigned short getNucleotide (unsigned long sequence, size_t position) {
  return ((sequence >> (position * NUCLEOTIDES_BITS)) & NUCLEOTIDE_MASK);
}

void printSequence (uint64_t sequence, size_t length) {
  for (unsigned int i = 0; i < length; ++i) {
    printf("%c", "ACGT"[sequence & NUCLEOTIDE_MASK]);
    sequence >>= NUCLEOTIDES_BITS;
  }
}

void printSequenceLong (unsigned short *sequence, size_t length) {
  for (unsigned int i = 0; i < length; ++i) {
    printf("%c", "ACGT"[sequence[i]]);
  }
}

#endif
