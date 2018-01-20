#!/usr/bin/env python2

import collections
import intervaltree
import re
import sys

# Configure the path to pybam
sys.path.append('./pybam/')
import pybam

gtf_file = sys.argv[1]
bam_files = sys.argv[2:]

gene_names = collections.defaultdict(intervaltree.IntervalTree)
exons = collections.defaultdict(intervaltree.IntervalTree)
with open(gtf_file, 'r') as f:
    for line in f:
        if "gene_type protein_coding" not in line:
            next
        t = line.strip().split("\t")
        if t[2] == "gene":
            i1 = t[8].find('gene_name')
            i2 = t[8].find(';', i1)
            gn = t[8][i1+9:i2].strip()
            gene_names[t[0]][int(t[3]):(int(t[4])+1)] = gn
        elif (t[2] == "exon") and ("transcript_type protein_coding" in line):
            exons[t[0]][int(t[3]):(int(t[4])+1)] = 1


bam_parsers = [pybam.read(bam_file, ['sam_qname', 'sam_rname', 'sam_pos1', 'sam_cigar_string', 'sam_cigar_list', 'sam_tags_list']) for bam_file in bam_files]


# Input: iterator (...data..., tags)
# Output: iterator (...data..., NM_tag)
# Description: Discards the read that are not uniquely mapped (i.e. don't
#              have NH:i:1). Also replaces the tag list with the value of
#              the NM tag.
def discard_non_unique_mappings(bam_parser):
    for read in bam_parser:
        NH_value = [t[2] for t in read[-1] if t[0] == 'NH'][0]  # Assumes the tag is always present
        if NH_value == 1:
            NM_value = [t[2] for t in read[-1] if t[0] == 'NM'][0]  # Assumes the tag is always present
            yield (read[:-1] + (NM_value,))

def only_exonic_mappings(bam_parser):
    for read in bam_parser:
        (read_name, chrom_name, start_pos, _, cigar_list, _) = read
        end_pos = start_pos + mapping_length(cigar_list) - 1
        if exons[chrom_name][start_pos:(end_pos+1)]:
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


def mapping_length(cigar_list):
    return sum(-n if c == 'I' else n for (n, c) in cigar_list)

def select_same_gene(group_iterator):
    for mappings in group_iterator:
        genes_seen = set()
        for pair in mappings[1]:
            if pair is not None:
                for (chrom_name, start_pos, _, cigar_list, _) in pair:
                    end_pos = start_pos + mapping_length(cigar_list) - 1
                    names_here = [i.data for i in gene_names[chrom_name][start_pos:(end_pos+1)]]
                    genes_seen.update(names_here)
        if len(genes_seen) == 1:
            yield mappings


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
for (read_name,mappings) in select_same_gene(merged_iterators([paired_reads_parser(discard_non_unique_mappings(p)) for p in bam_parsers])):
#for (names,(read_name,mappings)) in select_same_gene(merged_iterators([paired_reads_parser(discard_non_unique_mappings(p)) for p in bam_parsers])):
#for (read_name,mappings) in merged_iterators([paired_reads_parser(discard_non_unique_mappings(p)) for p in bam_parsers]):
    line = [read_name]
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

