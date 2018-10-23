// ==========================================================================
//                              srnaMapper
// ==========================================================================
// Copyright (C) 2018 Matthias Zytnicki, INRA
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
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
// Parse parameters
// ==========================================================================
#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <unordered_map>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

static const unsigned int N_NUCLEOTIDES           = 5;
static const unsigned int TRIPLET                 = 3;
static const unsigned int N_TRIPLETS              = std::pow(N_NUCLEOTIDES, TRIPLET);
static const unsigned int TRIPLET_MASK            = std::pow(N_NUCLEOTIDES, TRIPLET-1);
static const char         DNA5_TO_CHAR []         = { 'A', 'C', 'G', 'T', 'N' };

unsigned int getCode (char c) {
  switch (c) {
    case 'a':
    case 'A':
      return 0;
    case 'c':
    case 'C':
      return 1;
    case 'g':
    case 'G':
      return 2;
    case 't':
    case 'T':
    case 'u':
    case 'U':
      return 3;
    default:
      return 4;
  }
}

inline void updateQuality (std::string & s, std::string & t) {
    std::transform(s.begin(), s.end(), t.begin(), t.begin(), [] (char & c, char & d) { return std::max<char>(c, d); });
}

inline void updateQuality (seqan::CharString & s, seqan::CharString & t) {
    std::transform(begin(s), end(s), begin(t), begin(t), [] (char & c, char & d) { return std::max<char>(c, d); });
}

struct Parameters {
    std::vector < std::string >  readsFileNames;
    seqan::CharString outputFileName;
    unsigned int nReadsFiles;
    unsigned int lowComplexityThreshold;
    Parameters (): readsFileNames(), outputFileName(), nReadsFiles(0), lowComplexityThreshold(0) {}
};

Parameters parameters;

seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const ** argv) {
  seqan::ArgumentParser parser("srnaMapper");
  setShortDescription(parser, "Mapper for short short reads!");
  setVersion(parser, "0.0");
  seqan::addOption(parser, seqan::ArgParseOption("r", "reads", "The input reads file", seqan::ArgParseArgument::INPUT_FILE, "IN", true));
  seqan::addOption(parser, seqan::ArgParseOption("o", "output", "The output file", seqan::ArgParseArgument::OUTPUT_FILE, "IN"));
  seqan::addOption(parser, seqan::ArgParseOption("f", "filter", "Low complexity filter", seqan::ArgParseArgument::INTEGER, "INT"));
  setRequired(parser, "r");
  setRequired(parser, "o");
  setDefaultValue(parser, "filter", "8");
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK) {
    return res;
  }
  parameters.readsFileNames = getOptionValues(parser, "reads");
  getOptionValue(parameters.outputFileName, parser, "output");
  parameters.nReadsFiles = getOptionValueCount(parser, "reads");
  getOptionValue(parameters.lowComplexityThreshold, parser, "filter");
  return seqan::ArgumentParser::PARSE_OK;
}

class Tree {
private:
  size_t depth;
  std::vector < std::array < size_t, N_NUCLEOTIDES > > sequences;
  std::vector < std::string > qualities;
  std::vector < std::vector < unsigned int > > counts;
  std::unordered_map < size_t, size_t > sequence2quality;
  const static size_t NO_DATA = static_cast < size_t > (-1);
  std::vector < unsigned int > tripletCounts;

  void createCell () {
    sequences.emplace_back(std::array < size_t, N_NUCLEOTIDES > {NO_DATA, NO_DATA, NO_DATA, NO_DATA, NO_DATA});
  }

public:

  Tree (): depth(0), sequences(), qualities(), counts(), sequence2quality(), tripletCounts(N_TRIPLETS, 0) {
    createCell();
  }

  size_t getDepth () const {
    return depth;
  }

  void add (std::string & sequence, std::string & quality, size_t idReadsFile) {
    size_t       pos             = 0;
    size_t       size            = sequences.size();
    unsigned int tripletId       = 0;
    bool         isTriplet       = false;
    depth = std::max<size_t> (depth, sequence.size());
    std::fill(tripletCounts.begin(), tripletCounts.end(), 0);
    for (size_t i = sequence.size(); i > 0; --i) {
      char c = sequence[i-1];
      int  b = getCode(c);
      tripletId *= N_NUCLEOTIDES;
      tripletId += b;
      if ((! isTriplet) && (i <= sequence.size() - TRIPLET + 1)) {
        isTriplet = true;
      }
      if (isTriplet) {
        tripletCounts[tripletId] += 1;
        if (tripletCounts[tripletId] >= parameters.lowComplexityThreshold) {
          return;
        }
        tripletId %= TRIPLET_MASK;
      }
      if (sequences[pos][b] == NO_DATA) {
        createCell();
        pos = sequences[pos][b] = size;
        ++size;
      }
      else {
        pos = sequences[pos][b];
      }
    }
    auto it = sequence2quality.find(pos);
    if (it == sequence2quality.end()) {
      std::vector < unsigned int > count (parameters.nReadsFiles, 0);
      count[idReadsFile] = 1;
      qualities.push_back(quality);
      counts.push_back(count);
      sequence2quality[pos] = qualities.size() - 1;
    }
    else {
      updateQuality(qualities[(*it).second], quality);
      counts[(*it).second][idReadsFile] += 1;
    }
  }

  void print (unsigned int & readId, std::string & readString, size_t readPos, std::ostream & os, size_t treePos) const {
    auto   it      = sequence2quality.find(treePos);
    size_t nextPos = 0;
    if (it != sequence2quality.end()) {
      os << "@read_" << (++readId) << " x";
      for (unsigned int c: counts[(*it).second]) {
        os << "_" << c;
      }
      os << "\n" << (readString.data()+readPos+1) << "\n+\n" << qualities[(*it).second] << "\n";
    }
    for (size_t i = 0; i < N_NUCLEOTIDES; ++i) {
      nextPos = sequences[treePos][i];
      if (nextPos != NO_DATA) {
        readString[readPos] = DNA5_TO_CHAR[i];
        print(readId, readString, readPos-1, os, nextPos);
      }
    }
  }

  friend std::ostream & operator<<(std::ostream& os, const Tree & t);
};

std::ostream & operator<<(std::ostream & os, const Tree & t) {
  std::string readString (t.getDepth()+1, 0);
  unsigned int readId = 0;
  t.print(readId, readString, t.getDepth()-1, os, 0);
  return os;
}

int readReadsFile () {
  Tree tree;
  for (size_t idReadsFile = 0; idReadsFile < parameters.nReadsFiles; ++idReadsFile) {
    std::string &readsFileName = parameters.readsFileNames[idReadsFile];
    std::ifstream readsFile (seqan::toCString(readsFileName), std::ifstream::in);
    if (! readsFile.is_open()) {
      std::cerr << "ERROR: Could not open the file '" << readsFileName << "'.\n";
      return 1;
    }
    std::string line, sequence, quality;
    std::cerr << "Reading '" << readsFileName << "'...\n";
    while (std::getline(readsFile, line)) {
      std::getline(readsFile, sequence);
      std::getline(readsFile, line);
      std::getline(readsFile, quality);
      tree.add(sequence, quality, idReadsFile);
    }
    std::cerr << "... done.\n";
  }
  std::ofstream outputFile(seqan::toCString(parameters.outputFileName), std::ofstream::out);
  if (! outputFile.is_open()) {
    std::cerr << "ERROR: Could not open the file '" << parameters.outputFileName << "'.\n";
    return 1;
  }
  std::cerr << "Printing reads to '" << parameters.outputFileName << "'...\n";
  outputFile << tree;
  std::cerr << "... done.\n";
  return 0;
}

int main(int argc, char const ** argv) {
  seqan::ArgumentParser::ParseResult res = parseCommandLine(argc, argv);

  if (res != seqan::ArgumentParser::PARSE_OK) {
    return res == seqan::ArgumentParser::PARSE_ERROR;
  }
  int readRes = readReadsFile();
  if (readRes != 0) {
    return res;
  }
  return 0;
}
