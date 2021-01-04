#!/usr/bin/env python2

import collections
import gzip
import optparse
import os
import re
import sys
import time

import pybam


# Define a wrapper around quicksect.IntervalTree or intervaltree.IntervalTree, depending
# on which one is available, so that the rest of the script can work with both.
# quicksect is much faster than intervaltree, but harder to install.

# Throughout the script, start and end coordinates are all 1-based and both included
try:
    import quicksect

    class myintervaltree(quicksect.IntervalTree):

        def __init__(self):
            super(myintervaltree, self).__init__()
            self.interval_to_data = {}
            self.has_overlap = self.search

        # quicksect cannot associate intervals to data, so we need to it ourselves
        def add_data(self, start, end, data):
            interval = quicksect.Interval(start, end)
            self.insert( interval )
            self.interval_to_data[ interval ] = data

        def data_overlapping(self, start, end):
            return [self.interval_to_data[i] for i in self.search(start, end)]

except ImportError:
    import intervaltree

    class myintervaltree(intervaltree.IntervalTree):

        def __init__(self):
            super(myintervaltree, self).__init__()
            # Deal with versions 2.* and 3.*
            if hasattr(intervaltree.IntervalTree, "search"):
                self.search_func = self.search
            else:
                self.search_func = self.overlap
            self.has_overlap = self.search_func

        # In intervaltree.IntervalTree the intervals include the lower bound but not the upper bound.
        def add_data(self, start, end, data):
            self.addi(start, end+1, data)

        def has_overlap(self, start, end):
            return self.search_func(start, end+1)

        def data_overlapping(self, start, end):
            return [i.data for i in self.search_func(start, end+1)]


parser = optparse.OptionParser(usage = "usage: %prog [options] gtf_file bam_file_1 bam_file_2 ...")
parser.add_option("-r", "--rename_gene", nargs = 2, action = 'append', dest = 'gene_renames', metavar = 'NAME_TO_REPLACE NEW_NAME',
        help = 'Replace some erroneous gene names')
parser.add_option("-g", "--gene_types", dest = 'gene_types', default = 'protein_coding',
        help = 'Comma-separated list of gene biotypes to use [default: %default]. Use an empty string for no filtering')
parser.add_option("-t", "--transcript_types", dest = 'transcript_types', default = 'protein_coding',
        help = 'Comma-separated list of transcript biotypes to use for the exon-overlap filtering [default: %default]. Use an empty string for no filtering')
parser.add_option("-n", "--no_gtf_filter", action = 'store_true',
        help = 'Do not use a GTF file to filter the reads. The command line arguments are then expected to all be BAM files.')
(options, args) = parser.parse_args()

if options.no_gtf_filter:
    parser.usage = parser.usage.replace('gtf_file ', '')
    gtf_file = None
else:
    if len(args) == 0:
        parser.error("No GTF/BAM files given")
    gtf_file = args.pop(0)

if len(args) == 0:
    parser.error("No BAM files given")

bam_files = args
print >> sys.stderr, "BAM+GTF merger for AltHapAlign called with these options"
print >> sys.stderr, "\tgtf_file:", gtf_file
print >> sys.stderr, "\tbam_files:", bam_files
print >> sys.stderr, "\toptions:", options
print >> sys.stderr

ref_time = time.time()
def get_elapsed():
    time_now = time.time()
    global ref_time
    elapsed = time_now - ref_time
    ref_time = time_now
    return elapsed


if not options.no_gtf_filter:
    print >> sys.stderr, "Loading the GTF file ... ",
gene_type_filter = re.compile('gene_type "?(%s)"?;' % options.gene_types.replace(",", "|")) if options.gene_types else None
transcript_type_filter = re.compile('transcript_type "?(%s)"?;' % options.transcript_types.replace(",", "|")) if options.transcript_types else None
gene_renames = dict(options.gene_renames) if options.gene_renames else {}

gene_names = collections.defaultdict(myintervaltree)
exons = collections.defaultdict(myintervaltree)
n_genes = set()
n_exons = 0
if options.no_gtf_filter:
    gtf_file = os.devnull
with gzip.GzipFile(gtf_file, 'r') if gtf_file.endswith('.gz') else open(gtf_file, 'r') as f:
    for line in f:
        t = line.strip().split("\t")
        if gene_type_filter and not gene_type_filter.search(t[8]):
            continue
        if transcript_type_filter and not transcript_type_filter.search(t[8]):
            continue
        if t[2] == "exon":
            if 'gene_name' not in t[8]:
                print >> sys.stderr, "Failed\nERROR: no 'gene_name' attribute for this exon:", line,
                sys.exit(1)
            i1 = t[8].find('gene_name')
            i2 = t[8].find(';', i1)
            gn = t[8][i1+9:i2].strip().strip('"')
            name = gene_renames.get(gn, gn)
            start = int(t[3])
            end = int(t[4])
            gene_names[t[0]].add_data(start, end, name)
            n_genes.add(name)
            exons[t[0]].add_data(start, end, 1)
            n_exons += 1
if not options.no_gtf_filter:
    print >> sys.stderr, "Done (%.2f seconds): %d genes and %d exons" % (get_elapsed(), len(n_genes), n_exons)

print >> sys.stderr, "Opening the BAM files ...",
bam_parsers = [pybam.read(bam_file, ['sam_qname', 'sam_rname', 'sam_pos1', 'sam_cigar_list', 'sam_tags_list']) for bam_file in bam_files]
print >> sys.stderr, " Done (%.2f seconds)" % (get_elapsed(),)

n_bam_aligns = 0

# Input: BAM iterator (with tag list)
# Output: BAM iterator (with the value of the NM tag instead of the tag list)
# Description: Discards the read that are not uniquely mapped (i.e. don't
#              have NH:i:1). Also replaces the tag list with the value of
#              the NM tag.
def discard_non_unique_mappings(bam_parser):
    for read in bam_parser:
        global n_bam_aligns
        n_bam_aligns += 1
        NH_value = [t[2] for t in read[-1] if t[0] == 'NH'][0]  # Assumes the tag is always present
        if NH_value == 1:
            NM_value = [t[2] for t in read[-1] if t[0].upper() == 'NM'][0]  # Assumes the tag is always present
            yield (read[:-1] + (NM_value,))


# Input: BAM iterator
# Output: BAM iterator (with the end position (adjusted for the interval-tree searches) instead of the cigar alignment)
# Description: Discards the reads that do not overlap any exons
def only_exonic_mappings(bam_parser):
    for (read_name, chrom_name, start_pos, cigar_list, NM_value) in bam_parser:
        end_pos = start_pos + mapping_length(cigar_list) - 1
        if (chrom_name in exons) and exons[chrom_name].has_overlap(start_pos, end_pos):
            yield (read_name, chrom_name, start_pos, end_pos, NM_value)

# Input: BAM iterator
# Output: BAM iterator (with the end position (adjusted for the interval-tree searches) instead of the cigar alignment)
# Description: Equivalent of only_exonic_mappings that only transforms the data without doing any filtering
def all_mappings(bam_parser):
    for (read_name, chrom_name, start_pos, cigar_list, NM_value) in bam_parser:
        end_pos = start_pos + mapping_length(cigar_list) - 1
        yield (read_name, chrom_name, start_pos, end_pos, NM_value)


# Input: iterator (read_name, ...data...)
# Output: iterator (read_name, [(...data1...), (...data2...)])
# Description: Group consecutive pairs of reads that have the same name.
#              Singletons and reads mapping to 3 or more times are thus
#              discarded.  This assumes that the input BAM file is sorted
def paired_reads_parser(bam_parser):
    try:
        read = next(bam_parser)
        last_reads = [read[1:]]
        last_read_name = read[0]
        n_reads = 1
    except StopIteration:
        return
    for read in bam_parser:
        this_read_name = read[0]
        # Same read name as last time
        if last_read_name != this_read_name:
            if n_reads == 2:
                yield (last_read_name, last_reads)
            last_reads = []
            last_read_name = this_read_name
            n_reads = 0
        last_reads.append(read[1:])
        n_reads += 1


# Input: list of iterators (read_name, data)
# Output: iterator (read_name, [data1 | None, data2 | None, ...])
# Description: For each read name, groups the reads from each BAM file,
#              assuming that the files are sorted by read name
# Remarks: This is similar to a k-way merge algorithm. Although the merge
#          itself can be achieved in O(log(k)) in theory, my implementations
#          did not provide any improvements, especially because the output
#          of the function is O(k). So in the end I'm merging the streams in O(k)
def merged_iterators(input_parsers):
    current_reads = [parser.next() for parser in input_parsers]   # Assumes none of the iterators is empty
    active_bam_parsers = len(input_parsers)
    while active_bam_parsers:
        read_names = [cr[0] for cr in current_reads if cr[0] is not None]
        next_read_name = min(read_names)
        yield (next_read_name, [cr[1] if cr[0] == next_read_name else None for cr in current_reads])
        for (i, cr) in enumerate(current_reads):
            if cr[0] == next_read_name:
                try:
                    current_reads[i] = input_parsers[i].next()
                except StopIteration:
                    current_reads[i] = (None,)
                    active_bam_parsers -= 1


# Input: CIGAR alignment already parsed to a list of (number, character)
# Output: integer
# Description: Computes the length of the alignment on the genome
def mapping_length(cigar_list):
    s = 0
    for (n, c) in cigar_list:
        if c == 'I':
            s -= n
        else:
            s += n
    return s


# Input: iterator (read_name, [data1 | None, data2 | None, ...])
# Output: iterator (read_name, [data1 | None, data2 | None, ...])
# Description: For each group of reads, finds the gene names on the genome
#              (using the alignment) and only keeps the groups that map to
#              a single gene name.
def select_same_gene(group_iterator):
    for (read_name, mappings) in group_iterator:
        genes_seen = set()
        for pair in mappings:
            if pair is not None:
                for (chrom_name, start_pos, end_pos, NM_value) in pair:
                    names_here = gene_names[chrom_name].data_overlapping(start_pos, end_pos)
                    genes_seen.update(names_here)
        if len(genes_seen) == 1:
            yield (read_name, genes_seen.pop(), mappings)


def toString(group_iterator):
    for g in group_iterator:
        line = list(g[:-1])     # read_name, gene_name if a GTF file is given
        for m in g[-1]:         # mappings
            line.append('NA' if m is None else str(m[0][3] + m[1][3]))
        yield "\t".join(line)


n_groups = 0
last_n_bam_aligns = 0
partial_time = ref_time
print >> sys.stderr, "Reading the BAM files ..."
headers = ["read_name"]
if options.no_gtf_filter:
    it = toString(merged_iterators([paired_reads_parser(all_mappings(discard_non_unique_mappings(p))) for p in bam_parsers]))
else:
    headers = headers + ["gene_name"]
    it = toString(select_same_gene(merged_iterators([paired_reads_parser(only_exonic_mappings(discard_non_unique_mappings(p))) for p in bam_parsers])))
headers = headers + bam_files
print "\t".join(headers)
for s in it:
    print s
    n_groups += 1
    if not n_groups % 10000:
        print >> sys.stderr, "Found %d paired reads across all BAM files (%d raw reads processed -- %.2f per second)" % (n_groups, n_bam_aligns, (n_bam_aligns-last_n_bam_aligns)/get_elapsed())
        last_n_bam_aligns = n_bam_aligns
ref_time = partial_time
print >> sys.stderr, "Finished reading the BAM files in %.2f seconds: %d paired reads in total" % (get_elapsed(), n_groups)

