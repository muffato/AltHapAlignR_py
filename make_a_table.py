#!/usr/bin/env python2

import sys

# Configure the path to pybam
sys.path.append('./pybam/')
import pybam

bam_files = sys.argv[1:]

bam_parsers = [pybam.read(bam_file, ['sam_qname', 'sam_pos1', 'sam_cigar_string', 'sam_tags_list']) for bam_file in bam_files]


# Input: iterator (read_name, start_pos, cigar_string, tags)
# Output: iterator (read_name, start_pos, cigar_string, NM_tag)
# Description: Discards the read that are not uniquely mapped (i.e. don't
#              have NH:i:1). Also replaces the tag list with the value of
#              the NM tag.
def discard_non_unique_mappings(bam_parser):
    for read in bam_parser:
        NH_value = [t[2] for t in read[3] if t[0] == 'NH'][0]  # Assumes the tag is always present
        if NH_value == 1:
            NM_value = [t[2] for t in read[3] if t[0] == 'NM'][0]  # Assumes the tag is always present
            yield (read[:3] + (NM_value,))


# Input: iterator (read_name, ...data...)
# Output: iterator (read_name, (...data1...), (...data2...))
# Description: Group consecutive pairs of reads that have the same name.
#              Singletons and reads mapping to 3 or more times are thus
#              discarded.  This assumes that the input BAM file is sorted
def paired_reads_parser(bam_parser):
    last_read = None
    read_count = 0
    for read in bam_parser:
        # Same read as last time
        if (last_read is None) or (last_read[0] != read[0]):
            if read_count == 2:
                yield (read[0], last_read[1:], read[1:])
            read_count = 0
        last_read = read
        read_count += 1


# Input: list of iterators (read_name, (...data1...), (...data2...))
# Output: iterator (read_name, [iterator_index, (...data1...), (...data2...)})
# Description: For each read name, groups the reads from each BAM file,
#              assuming that the files are sorted by read name
def merged_iterators(input_parsers):
    current_reads = [parser.next() for parser in input_parsers]   # Assumes none of the iterators is empty
    active_bam_parsers = len(input_parsers)
    while active_bam_parsers:
        read_names = [cr[0] for cr in current_reads if cr[0] is not None]
        next_read_name = min(read_names)
        yield (next_read_name, [cr[1:] if cr[0] == next_read_name else None for cr in current_reads])
        for (i, cr) in enumerate(current_reads):
            if cr[0] == next_read_name:
                try:
                    current_reads[i] = input_parsers[i].next()
                except StopIteration:
                    current_reads[i] = (None,)
                    active_bam_parsers -= 1


n = 0
#for (read_name, start_pos, cigar_string, tags) in bam_parsers[0]:
    #print (read_name, start_pos, cigar_string, tags)
#for g in paired_reads_parser(discard_non_unique_mappings(bam_parsers[0])):
#for g in merged_iterators([paired_reads_parser(discard_non_unique_mappings(p)) for p in bam_parsers]):
    #print g
    #n += 1
    #if n == 5:
        #break

for (read_name,mappings) in merged_iterators([paired_reads_parser(discard_non_unique_mappings(p)) for p in bam_parsers]):
    line = [read_name]
    for m in mappings:
        if m is None:
            line.extend(['NA'] * 6)
        else:
            line.extend(m[0])
            line.extend(m[1])
    print "\t".join(map(str,line))
    #n += 1
    #if n == 5:
        #break

