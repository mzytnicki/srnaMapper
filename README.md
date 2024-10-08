# Install

## Manual install

* Download the tool with `git clone --recurse-submodules https://github.com/mzytnicki/srnaMapper.git`
* Go the directory with `cd srnaMapper`
* Compile with `make`

You should have the `zlib` library, and a C++11 compiler.

## With bioconda

If you have [bioconda](https://bioconda.github.io/) ready, you can install srnaMapper with:

    conda install -c bioconda srnamapper


# Generating the genome index

srnaMapper uses the `bwa` API, and the `bwa` index files.
Use `bwa index` to generate the index files.

## Note

As in `bwa`, ambiguous nucleotides in the genome are converted to a random nucleotide.

The procedure is identical for the reads. 

# Using srnaMapper

Once compiled type `./srnaMapper` *parameters*.

Compulsory parameters:

* `-r` *string*: file name in FASTQ format
* `-g` *string*: prefix of the genome database (produced by `bwa build`)
* `-o` *string*: output file in SAM format

Optional parameters:

* `-e` *int*: maximum number of errors (default: 2)
* `-t` *int*: number of threads (default: 1)
* `-n` *int*: discard reads when they map more than n times (default: 5)
* `-f` *int*: low complexity threshold, more is more lenient (default: 6)
* `-u`: if set, print all the mapped reads in a unique SAM file (with the counts for each sample)
* `-s` *int*: set the random seed (time otherwise)
* `-h`: the help message

Notes:

* The `-r` option should be repeated once per input file.
* Unless the the `-u` option is set, the `-o` option should also be repeated once per input file.

Example:

    ./srnaMapper -r cond1_rep1.fastq -r cond1_rep2.fastq -r cond2_rep1.fastq -r cond2_rep2.fastq -g genome -o cond1_rep1.sam -o cond1_rep2.sam -o cond2_rep1.sam -o cond2_rep2.sam -t 10
