#!/usr/bin/env python2

import collections
import intervaltree
import optparse
import re
import sys

# quicksect is much faster than intervaltree, but harder to install
has_quicksect = True
try:
    import quicksect
except ImportError:
    has_quicksect = False

# Configure the path to pybam
sys.path.append('./pybam/')
import pybam

parser = optparse.OptionParser(usage = "usage: %prog [options] gtf_file bam_file_1 bam_file_2 ...")
parser.add_option("-r", "--rename_gene", nargs = 2, action = 'append', dest = 'gene_renames', metavar = 'NAME_TO_REPLACE NEW_NAME',
        help = 'Replace some erroneous gene names')
parser.add_option("-g", "--gene_types", dest = 'gene_types', default = 'protein_coding',
        help = 'Comma-separated list of gene biotypes to use [default: %default]. Use an empty string for no filtering')
parser.add_option("-t", "--transcript_types", dest = 'transcript_types', default = 'protein_coding',
        help = 'Comma-separated list of transcript biotypes to use for the exon-overlap filtering [default: %default]. Use an empty string for no filtering')
(options, args) = parser.parse_args()

gtf_file = args[0]
bam_files = args[1:]


gene_type_filter = re.compile("gene_type (%s);" % options.gene_types.replace(",", "|")) if options.gene_types else None
transcript_type_filter = re.compile("transcript_type (%s);" % options.transcript_types.replace(",", "|")) if options.transcript_types else None
gene_renames = dict(options.gene_renames) if options.gene_renames else {}

gene_names = collections.defaultdict(intervaltree.IntervalTree)
exons = collections.defaultdict(quicksect.IntervalTree if has_quicksect else intervaltree.IntervalTree)
with open(gtf_file, 'r') as f:
    for line in f:
        t = line.strip().split("\t")
        if gene_type_filter and not gene_type_filter.search(t[8]):
            continue
        if t[2] == "gene":
            i1 = t[8].find('gene_name')
            i2 = t[8].find(';', i1)
            gn = t[8][i1+9:i2].strip()
            gene_names[t[0]][int(t[3]):(int(t[4])+1)] = gene_renames.get(gn, gn)
        elif t[2] == "exon":
            if transcript_type_filter and not transcript_type_filter.search(t[8]):
                continue
            if has_quicksect:
                exons[t[0]].add(int(t[3]), int(t[4])+1)
            else:
                exons[t[0]][int(t[3]):(int(t[4])+1)] = 1


bam_parsers = [pybam.read(bam_file, ['sam_qname', 'sam_rname', 'sam_pos1', 'sam_cigar_string', 'sam_cigar_list', 'sam_tags_list']) for bam_file in bam_files]


# Input: BAM iterator (with tag list)
# Output: BAM iterator (with the value of the NM tag instead of the tag list)
# Description: Discards the read that are not uniquely mapped (i.e. don't
#              have NH:i:1). Also replaces the tag list with the value of
#              the NM tag.
def discard_non_unique_mappings(bam_parser):
    for read in bam_parser:
        NH_value = [t[2] for t in read[-1] if t[0] == 'NH'][0]  # Assumes the tag is always present
        if NH_value == 1:
            NM_value = [t[2] for t in read[-1] if t[0] == 'NM'][0]  # Assumes the tag is always present
            yield (read[:-1] + (NM_value,))


# Input: BAM iterator
# Output: BAM iterator
# Description: Discards the reads that do not overlap any exons
def only_exonic_mappings(bam_parser):
    for read in bam_parser:
        (read_name, chrom_name, start_pos, _, cigar_list, _) = read
        end_pos_plus_1 = start_pos + mapping_length(cigar_list)
        if exons[chrom_name].search(start_pos, end_pos_plus_1):
            yield read


# Input: iterator (read_name, ...data...)
# Output: iterator (read_name, [(...data1...), (...data2...)])
# Description: Group consecutive pairs of reads that have the same name.
#              Singletons and reads mapping to 3 or more times are thus
#              discarded.  This assumes that the input BAM file is sorted
def paired_reads_parser(bam_parser):
    last_reads = []
    last_read_name = None
    for read in bam_parser:
        # Same read as last time
        if (len(last_reads) == 0) or (last_read_name != read[0]):
            if len(last_reads) == 2:
                yield (last_read_name, last_reads)
            last_reads = []
            last_read_name = read[0]
        last_reads.append(read[1:])


# Input: list of iterators (read_name, data)
# Output: iterator (read_name, [data1 | None, data2 | None, ...])
# Description: For each read name, groups the reads from each BAM file,
#              assuming that the files are sorted by read name
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
    return sum(-n if c == 'I' else n for (n, c) in cigar_list)


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
                for (chrom_name, start_pos, _, cigar_list, _) in pair:
                    end_pos_plus_1 = start_pos + mapping_length(cigar_list)
                    names_here = [i.data for i in gene_names[chrom_name][start_pos:end_pos_plus_1]]
                    genes_seen.update(names_here)
        if len(genes_seen) == 1:
            yield (read_name, genes_seen.pop(), mappings)


n = 0
#for (read_name, start_pos, cigar_string, tags) in bam_parsers[0]:
    #print (read_name, start_pos, cigar_string, tags)
#for g in paired_reads_parser(discard_non_unique_mappings(bam_parsers[0])):
#for g in merged_iterators([paired_reads_parser(discard_non_unique_mappings(p)) for p in bam_parsers]):
    #print g
    #n += 1
    #if n == 5:
        #break

#for (read_name,mappings) in select_same_gene(merged_iterators([paired_reads_parser(only_exonic_mappings(discard_non_unique_mappings(p))) for p in bam_parsers])):
for (read_name, gene_name, mappings) in select_same_gene(merged_iterators([paired_reads_parser(only_exonic_mappings(discard_non_unique_mappings(p))) for p in bam_parsers])):
#for (names,(read_name,mappings)) in select_same_gene(merged_iterators([paired_reads_parser(discard_non_unique_mappings(p)) for p in bam_parsers])):
#for (read_name,mappings) in merged_iterators([paired_reads_parser(discard_non_unique_mappings(p)) for p in bam_parsers]):
    line = [read_name, gene_name]
    for m in mappings:
        if m is None:
            line.extend(['NA'] * 6)
        else:
            line.extend(m[0][x] for x in (1,2,4))
            line.extend(m[1][x] for x in (1,2,4))
    #print names
    print "\t".join(map(str,line))
    #print
    #n += 1
    #if n == 100000:
        #break

