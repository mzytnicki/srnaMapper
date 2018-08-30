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
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <unordered_map>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#define NUMBER_OF_READ_PER_BATCH 100000
#define N_NUCLEOTIDES seqan::ValueSize<seqan::Dna5>::VALUE

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
    Parameters (): readsFileNames(), outputFileName(), nReadsFiles(0) {}
};

Parameters parameters;

seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const ** argv) {
  seqan::ArgumentParser parser("srnaMapper");
  setShortDescription(parser, "Mapper for short short reads!");
  setVersion(parser, "0.0");
  seqan::addOption(parser, seqan::ArgParseOption("r", "reads", "The input reads file", seqan::ArgParseArgument::INPUT_FILE, "IN", true));
  seqan::addOption(parser, seqan::ArgParseOption("o", "output", "The output file", seqan::ArgParseArgument::OUTPUT_FILE, "IN"));
  setRequired(parser, "r");
  setRequired(parser, "o");
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK) {
    return res;
  }
  parameters.readsFileNames = getOptionValues(parser, "reads");
  getOptionValue(parameters.outputFileName, parser, "output");
  parameters.nReadsFiles = getOptionValueCount(parser, "reads");
  return seqan::ArgumentParser::PARSE_OK;
}

class Tree {
private:
  std::vector < std::array < size_t, N_NUCLEOTIDES > > sequences;
  std::vector < std::string > qualities;
  std::vector < std::vector < unsigned int > > counts;
  std::unordered_map < size_t, size_t > sequence2quality;
  const static size_t NO_DATA = static_cast < size_t > (-1);

  void createCell () {
    sequences.emplace_back(std::array < size_t, N_NUCLEOTIDES > {NO_DATA, NO_DATA, NO_DATA, NO_DATA, NO_DATA});
  }

public:

  Tree (): sequences(), qualities(), counts(), sequence2quality() {
    createCell();
  }

  void add (std::string & sequence, std::string & quality, size_t idReadsFile) {
    size_t pos = 0;
    size_t size = sequences.size();
    for (char c: sequence) {
      int b = getCode(c);
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

  void print (unsigned int & count, std::string & currentString, size_t currentPos, std::ostream & os, size_t pos) const {
    auto it = sequence2quality.find(pos);
    if (it != sequence2quality.end()) {
      currentString[currentPos] = '\0';
      os << "@read_" << (++count) << " x";
      for (unsigned int c: counts[(*it).second]) {
        os << "_" << c;
      }
      os << "\n" << currentString.data() << "\n+\n" << qualities[(*it).second] << "\n";
    }
    for (size_t i = 0; i < N_NUCLEOTIDES; ++i) {
      size_t nextPos = sequences[pos][i];
      if (nextPos != NO_DATA) {
        char currentChar = seqan::Dna5(i);
        if (currentString.size() == currentPos) {
          currentString += currentChar;
        }
        else {
          currentString[currentPos] = currentChar;
        }
        print(count, currentString, currentPos+1, os, nextPos);
      }
    }
  }

  friend std::ostream & operator<<(std::ostream& os, const Tree & t);
};

std::ostream & operator<<(std::ostream & os, const Tree & t) {
  std::string currentString;
  unsigned int count = 0;
  t.print(count, currentString, 0, os, 0);
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
  outputFile << tree;
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
