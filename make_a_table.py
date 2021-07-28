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
            n_exons += 1
if not options.no_gtf_filter:
    print >> sys.stderr, "Done (%.2f seconds): %d genes and %d exons" % (get_elapsed(), len(n_genes), n_exons)

print >> sys.stderr, "Opening the BAM files ...",
bam_parsers = [pybam.read(bam_file, ['sam_qname', 'sam_rname', 'sam_pos1', 'sam_cigar_list', 'sam_tags_list']) for bam_file in bam_files]
print >> sys.stderr, " Done (%.2f seconds)" % (get_elapsed(),)


# Input: BAM read with tag list, and tag name
# Output: the value of the tag
# Exception: exit if the tag is not present, or in multiple copies
def get_tag_value(read, tag):
    values = [t[2] for t in read[-1] if t[0].upper() == tag]
    if not values:
        print >> sys.stderr, "ERROR: Could not find any %s tag in" % tag, read
        sys.exit(1)
    if len(values) > 1:
        print >> sys.stderr, "ERROR: More than 1 %s tag in" % tag, read
        sys.exit(1)
    return values[0]


# Input: BAM iterator with tag list
# Output: BAM iterator with selected tag values (NH and NM)
# Description: Extract the NH and NM values from the BAM alignments
n_bam_aligns = 0
def extract_tags(bam_parser):
    for read in bam_parser:
        global n_bam_aligns
        n_bam_aligns += 1
        NH_value = get_tag_value(read, 'NH')
        NM_value = get_tag_value(read, 'NM')
        yield (read[:-1] + (NH_value, NM_value))


# Input: iterator over paired alignments
# Output: iterator over a subset of the input
# Description: Discards the alignments that are not uniquely mapped pairs,
#              i.e. are singletons or have multiple hits / don't have NH:i:1.
n_singletons = 0
n_paired_alignments = 0
n_multiple_hits = 0
def select_paired_alignments(alignment_iterator):
    global n_singletons
    global n_paired_alignments
    global n_multiple_hits
    for data in alignment_iterator:
        n_alignments = len(data[1])
        NH_values = [alignment[-2] for alignment in data[1]]
        if n_alignments == 1:
            assert NH_values[0] == 1, data[0]
            n_singletons += 1
        elif n_alignments == 2:
            # Check that NH = 1 for both
            if NH_values[0] != 1:
                assert NH_values[1] != 1
                # Could be classified as singletons too
                n_multiple_hits += 2
            else:
                assert NH_values[1] == 1
                n_paired_alignments += 1
                yield data
        else:
            for NH_value in NH_values:
                assert NH_value != 1
            n_multiple_hits += n_alignments


# Input: iterator over paired alignments
# Output: iterator that has gene names too, but on a subset of the alignments
# Description: Discards the reads that do not overlap any exons
n_non_exonic_bam_aligns = 0
n_different_genes_pair = 0
n_partially_exonic = 0
def only_exonic_mappings(alignment_iterator):
    for (read_name, alignments) in alignment_iterator:
        ok = 0
        genes_seen = set()
        for (chrom_name, start_pos, cigar_list, NH_value, NM_value) in alignments:
            if chrom_name in gene_names:
                end_pos = start_pos + mapping_length(cigar_list) - 1
                names_here = gene_names[chrom_name].data_overlapping(start_pos, end_pos)
                if names_here:
                    genes_seen.update(names_here)
                    ok += 1
        if ok == 0:
            global n_non_exonic_bam_aligns
            n_non_exonic_bam_aligns += 1
        elif len(genes_seen) > 1:
            global n_different_genes_pair
            n_different_genes_pair += 1
        else:
            if ok == 1:
                global n_partially_exonic
                n_partially_exonic += 1
            yield (read_name, (genes_seen.pop(), alignments))


# Input: iterator (read_name, ...data...)
# Output: iterator (read_name, [(...data1...), (...data2...)])
# Description: Group consecutive pairs of reads that have the same name.
#              This assumes that the input BAM file is sorted.
def group_read_alignments(bam_parser):
    try:
        read = next(bam_parser)
        last_reads = [read[1:]]
        last_read_name = read[0]
    except StopIteration:
        return
    for read in bam_parser:
        this_read_name = read[0]
        # Same read name as last time
        if last_read_name != this_read_name:
            yield (last_read_name, last_reads)
            last_reads = []
            last_read_name = this_read_name
        last_reads.append(read[1:])
    yield (last_read_name, last_reads)


# Input: list of iterators (read_name, data)
# Output: iterator (read_name, [data1 | None, data2 | None, ...])
# Description: For each read name, groups the reads from each BAM file,
#              assuming that the files are sorted by read name
# Remarks: This is similar to a k-way merge algorithm. Although the merge
#          itself can be achieved in O(log(k)) in theory, my implementations
#          did not provide any improvements, especially because the output
#          of the function is O(k). So in the end I'm merging the streams in O(k)
n_groups = 0
def merged_iterators(input_parsers):
    current_reads = [parser.next() for parser in input_parsers]   # Assumes none of the iterators is empty
    active_bam_parsers = len(input_parsers)
    while active_bam_parsers:
        read_names = [cr[0] for cr in current_reads if cr[0] is not None]
        next_read_name = min(read_names)
        global n_groups
        n_groups += 1
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
n_unique_groups = 0
n_best_groups = 0
n_ambiguous_groups = 0
def select_same_gene(group_iterator):
    for (read_name, mappings) in group_iterator:
        genes_seen = {}
        for pair in mappings:
            if pair is not None:
                gene_name = pair[0]
                NM_score = sum(alignment[-1] for alignment in pair[1])
                if gene_name in genes_seen:
                    if NM_score < genes_seen[gene_name]:
                        genes_seen[gene_name] = NM_score
                else:
                    genes_seen[gene_name] = NM_score
        if len(genes_seen) == 1:
            global n_unique_groups
            n_unique_groups += 1
            yield (read_name, 'unique', genes_seen.keys()[0], mappings)
        else:
            best_NM = min(genes_seen.values())
            best_genes = [gene_name for (gene_name, NM_score) in genes_seen.items() if NM_score == best_NM]
            if len(best_genes) == 1:
                global n_best_groups
                n_best_groups += 1
                status = 'best'
            else:
                global n_ambiguous_groups
                n_ambiguous_groups += 1
                status = 'ambiguous'
            for gene_name in best_genes:
                yield (read_name, status, gene_name, [pair if pair and pair[0] == gene_name else None for pair in mappings])


def toString(g):
    line = list(data[:-1])     # read_name (always), gene_name (optional)
    for m in data[-1]:         # mappings
        line.append('NA' if m is None else str(m[1][0][-1] + m[1][1][-1]))
    return "\t".join(line)


last_n_bam_aligns = 0
partial_time = ref_time
print >> sys.stderr, "Reading the BAM files ..."
headers = ["read_name"]
bam_filters = [select_paired_alignments(group_read_alignments(extract_tags(bp))) for bp in bam_parsers]
if options.no_gtf_filter:
    it = merged_iterators(bam_filters)
else:
    headers = headers + ["gene_name"]
    it = select_same_gene(merged_iterators([only_exonic_mappings(bf) for bf in bam_filters]))
headers = headers + bam_files
print "\t".join(headers)
last_n_groups = 0
n_lines = 0
for data in it:
    print toString(data)
    n_lines += 1
    if n_groups >= last_n_groups+10000:
        print >> sys.stderr, "Found %d reads across all BAM files (%d alignments processed -- %.2f per second)" % (n_groups, n_bam_aligns, (n_bam_aligns-last_n_bam_aligns)/get_elapsed())
        last_n_bam_aligns = n_bam_aligns
        last_n_groups = n_groups

ref_time = partial_time
print >> sys.stderr, "Finished reading the BAM files in %.2f seconds" % get_elapsed()
print >> sys.stderr, "%d alignments across all %d BAM files" % (n_bam_aligns, len(bam_files))
print >> sys.stderr, "\t%d discarded (%.2f%%) - singletons " % (n_singletons, 100.*n_singletons/n_bam_aligns)
print >> sys.stderr, "\t%d discarded (%.2f%%) - multiple hits (NH != 1)" % (n_multiple_hits, 100.*n_multiple_hits/n_bam_aligns)
print >> sys.stderr, "%d paired alignments" % n_paired_alignments
print >> sys.stderr, "\t%d discarded (%.2f%%) - both not exonic" % (n_non_exonic_bam_aligns, 100.*n_non_exonic_bam_aligns/n_paired_alignments)
print >> sys.stderr, "\t%d discarded (%.2f%%) - in different genes" % (n_different_genes_pair, 100.*n_different_genes_pair/n_paired_alignments)
print >> sys.stderr, "\t%d kept (%.2f%%) - only one not exonic" % (n_partially_exonic, 100.*n_partially_exonic/n_paired_alignments)
print >> sys.stderr, "%d reads after grouping the BAM Files" % n_groups
if n_groups and (not options.no_gtf_filter):
    print >> sys.stderr, "Gene name assignment statistics"
    print >> sys.stderr, "\t%d reads (%.2f%%): all BAM files agree" % (n_unique_groups, 100.*n_unique_groups/n_groups)
    if n_best_groups:
        print >> sys.stderr, "\t%d reads (%.2f%%): multiple candidates, lowest NM score selected" % (n_best_groups, 100.*n_best_groups/n_groups)
    if n_ambiguous_groups:
        print >> sys.stderr, "\t%d reads (%.2f%%): multiple candidates, tie - %d candidates listed (%.2f per read on average)" % (n_ambiguous_groups, 100.*n_ambiguous_groups/n_groups, n_lines-n_unique_groups-n_best_groups, float(n_lines-n_unique_groups-n_best_groups)/n_ambiguous_groups)

#Return a non-zero code if we couldn't find any groups
if not n_groups:
    sys.exit(1)
