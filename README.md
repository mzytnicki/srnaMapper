# Compiling

Simply type `make`.


# Generating the genome index

srnaMapper uses the `bwa` API, and the `bwa` index files.
Use `bwa index` to generate the index files.


# Using srnaMapper

Once compiled type `./srnaMapper` *options*.

Compulsory options are:

* `-r` *reads*: the reads file, in FASTQ format
* `-g` *genome*: the genome index prefix (see previous section)
* `-o` *output*: the output file, in SAM format

Optional (?) options are:

* `-c` *filename*: the collapsed reads
* `-f` *filter*: the filter strenght used
* `-e` *n_errors*: the maximum number of errors
* `-h`: the help message
