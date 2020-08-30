#!/usr/bin/env python3
#
# splits up the name-sorted bamfile into one bam per contig
# does this in a single pass, opening as many files as there are contigs
# may be necessary to set ulimit -n higher on some systems
#
import pysam
import sys
import gfapy
import os
import json

gfa = gfapy.Gfa.from_file(sys.argv[1])
bam = pysam.AlignmentFile(sys.argv[2], "rb" )
contig_map = sys.argv[3]

out_bams = {}
out_read_count = {}
id = 0
contig_IDs = {}
for seg in gfa.segments:
    seg_bam = 'contig_'+str(id)+'.bam'
    contig_IDs[seg.name]='contig_'+str(id)
    id+=1
    out_bams[seg.name] = pysam.AlignmentFile(seg_bam, "wb", template=bam)
    out_read_count[seg.name] = 0

for read in bam.fetch(until_eof=True):
    if read.reference_name is None: continue
    out_bams[read.reference_name].write(read)
    out_read_count[read.reference_name] += 1

for contig in out_bams:
    out_bams[contig].close()
    if out_read_count[contig] == 0:
        del_bam = contig_IDs[contig]+'.bam'
        os.remove(del_bam)

with open(contig_map, 'w') as json_file:
    json.dump(contig_IDs, json_file)
