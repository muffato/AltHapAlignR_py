
## Installation

The python script has a few dependencies:

* [pybam](https://github.com/JohnLonginotto/pybam): "Pure Python" -but
  fast- library to read BAM files
* [intervaltree](https://pypi.python.org/pypi/intervaltree): "Pure Python"
  library that implements [interval trees](https://en.wikipedia.org/wiki/Interval_tree)
* [quicksect](https://pypi.python.org/pypi/quicksect): C/Python library
  that implements [interval trees](https://en.wikipedia.org/wiki/Interval_tree)
  too but is about 4x faster than `intervaltree`. Note that its
  installation may require [Cython](https://pypi.python.org/pypi/Cython)
  and a compiler (e.g. gcc) setup.

Only one of the last two is needed, `quicksect` being the preferred
option for performance reasons.

There are several ways of bringing them in, the easiest being with `pip`.
Note that you may want to first setup a [virtualenv](https://virtualenv.pypa.io)
before installing the dependencies, to ensure your environment is clean and
self-contained. For instance:

```sh
# Where the files are going to be stored
ALTHAPALIGN_VENV=$PWD/althapalign_virtualenv

# To create a "virtualenv" (only the first time)
virtualenv $ALTHAPALIGN_VENV

# To start using the "virtualenv"
source $ALTHAPALIGN_VENV/bin/activate

# To install the dependencies
pip install https://github.com/muffato/pybam/zipball/master
pip install cython
pip install quicksect

# To stop using it, once finished
deactivate
```


## Execution

### Inputs

* GTF file that covers all the haplotypes
* BAM files (1 per region). They *must* be sorted by *query name* with
  [Picard](http://broadinstitute.github.io/picard/). *Do not use samtools*
  to sort the files as it [does not follow the lexicographic
  order](https://github.com/samtools/hts-specs/issues/5).

### Command line

The general syntax is

```
./make_a_table.py [-r <wrong_gene_name> <correct_gene_name>] <path/to/GTF_file[.gz]> <path/to/BAM_file_1> <path/to/BAM_file_2> ...
```
The `-r` option can be repeated if several genes have the wrong name. For
instance, we use [Gencode
21](https://www.gencodegenes.org/releases/21.html) in the paper, which has
to be _fixed_ with `-r VARSL VARS2 -r C6orf205 MUC21`

### Output

The output is a tab-separated file with _n_+2 columns, _n_ being the number
of BAM files.

| read name | gene name | edit distance in BAM file 1 | edit distance in BAM file 2 | ... |
| --------- | --------- | --------------------------- | --------------------------- | --- |

Each edit distance is actually the _sum_ of the edit distances of the two
reads in each pair.

### Runtime statistics

The script's CPU-time is linear in the total number of reads
found in the BAM files. Depending on the CPU, it can parse between 30,000
and 80,000 reads per second. The analysis from the paper over 8 BAM files
comprising more than 13 million reads takes less than 3 minutes on an Intel i7-7500U
CPU @ 2.70GHz (on 1 core only). The memory consumption is limited than 20MB
with `quicksect`, 35MB with `intervaltree`, regardless of the number of reads.
It only depends on the number of genes and exons in the GTF file.

### Example command-line

```
./make_a_table.py -r VARSL VARS2 -r C6orf205 MUC21 gencode.v21.only_MHC.annotation.gtf.gz */*picard*bam > output_file
```
