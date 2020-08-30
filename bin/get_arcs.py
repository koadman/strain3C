#!/usr/bin/env python3
import pysam
import sys
import gfapy
import networkx as nx
from pysam import VariantFile
import json
import numpy

CONTIG_MARGIN = 200
MIN_SINGLECOPY_LENGTH = 20000
target_left = 0
target_right = 500000000

gfa = gfapy.Gfa.from_file(sys.argv[1])
vcf_in = pysam.VariantFile(sys.argv[2])
bam = pysam.AlignmentFile(sys.argv[3], "rb" )
id_map_file = open(sys.argv[4], 'r')
target_contig_id = sys.argv[5]

# find the target contig
id_map_json = json.load(id_map_file)
for contig in id_map_json:
    if id_map_json[contig] == target_contig_id:
        target_contig = contig

# get the reference sequences
ref_seqs = {}
for seg in gfa.segments:
    ref_seqs[seg.name] = seg.sequence

# get the variant sites and assign numeric IDs for use as array indices later
variants = {}
variant_IDs = {}
vID_lookup = []
vid = 0
for rec in vcf_in.fetch():
    if target_contig is not None and rec.contig != target_contig: continue
#    if rec.pos > 300: continue
    if rec.pos < target_left or rec.pos > target_right: continue
    variants[(rec.contig,rec.pos)] = [rec.ref, rec.alts]
    variant_IDs[(rec.contig,rec.pos)] = vid
    vID_lookup.append((rec.contig,rec.pos))
    vid += 1

print("found "+str(len(variants))+" variants")

local_var_counts = {}
hic_var_counts = {}
HIC_MIN_DIST = 800 # reads separated by more than this many bp are treated as Hi-C pairs
read_iter = bam.fetch(until_eof=True)
next_read = next(read_iter)
hic_pairs = 0
local_pairs = 0
single_reads = 0
try:
#if True:
    while True:
        # get a pair of aligned reads
        cur_read_1 = next_read
        next_read = next(read_iter)
        if cur_read_1.reference_name is None: continue
        if target_contig is not None and cur_read_1.reference_name != target_contig: continue
        if next_read.query_name == cur_read_1.query_name:
            cur_read_2 = next_read
            next_read = next(read_iter)
            if abs(cur_read_2.reference_start - cur_read_1.reference_start) > HIC_MIN_DIST:
                hic_pairs+=1
            else:
                local_pairs+=1
        else:
            cur_read_2 = None
            single_reads += 1

        # identify all the variant sites covered by the aligned reads
        var_sites = []
        for pair in cur_read_1.get_aligned_pairs(matches_only=True):
            if (cur_read_1.reference_name,pair[1]+1) in variants:
                refalt = 0
                if ref_seqs[cur_read_1.reference_name][pair[1]] != cur_read_1.query_sequence[pair[0]]:
                    refalt = 1
                var_sites.append((cur_read_1.reference_name,pair[1],refalt))
        if cur_read_2 is not None:
            for pair in cur_read_2.get_aligned_pairs(matches_only=True):
                if (cur_read_2.reference_name,pair[1]+1) in variants:
                    refalt = 0
                    if ref_seqs[cur_read_2.reference_name][pair[1]] != cur_read_2.query_sequence[pair[0]]:
                        refalt = 1
                    var_sites.append((cur_read_2.reference_name,pair[1],refalt))

        # tally up a co-observation for each pair of variant sites
        for p_i in range(len(var_sites)):
            for p_j in range(p_i+1,len(var_sites)):
                var_pair_ij = (var_sites[p_i][0],var_sites[p_i][1],var_sites[p_j][0],var_sites[p_j][1])
                indy = var_sites[p_i][2] + 2*var_sites[p_j][2]
                if var_pair_ij[0] == var_pair_ij[2] and abs(var_pair_ij[3]-var_pair_ij[1]) < HIC_MIN_DIST:
                    if not var_pair_ij in local_var_counts:
                        local_var_counts[var_pair_ij] = [0,0,0,0]
                    local_var_counts[var_pair_ij][indy] += 1
                else:
                    if not var_pair_ij in hic_var_counts:
                        hic_var_counts[var_pair_ij] = [0,0,0,0]
                    hic_var_counts[var_pair_ij][indy] += 1
except:
    pass
print("found "+str(len(local_var_counts))+" local variant site pairs in mapped reads")
print("found "+str(len(hic_var_counts))+" long range (HiC) variant site pairs in mapped reads")

print("hic pairs: "+str(hic_pairs))
print("local pairs: "+str(local_pairs))
print("single reads (no pair on target contig): "+str(single_reads))
dat = {
    'K': 0,
    'V': vid,
    'L_hic': len(hic_var_counts),
    'L_local': len(local_var_counts),
    'hic_linkcounts': [[0] * 4 for i in range(len(hic_var_counts))],
    'hic_linksites': [[0] * 2 for i in range(len(hic_var_counts))],
    'local_linkcounts': [[0] * 4 for i in range(len(local_var_counts))],
    'local_linksites': [[0] * 2 for i in range(len(local_var_counts))],
    'site_map': [[0] * 2 for i in range(len(vID_lookup))],
    'abundance_prior': 1000,
    'subsample': 10
}

l = 0
hic_minor = 0
for var in hic_var_counts:
    for i in range(4): dat['hic_linkcounts'][l][i] = hic_var_counts[var][i]
    gt1 = 0
    for i in range(1,4):
        if dat['hic_linkcounts'][l][i] > 0: gt1+=1
    if gt1 > 0:
        hic_minor += 1
    dat['hic_linksites'][l][0] = variant_IDs[(var[0],var[1]+1)]+1
    dat['hic_linksites'][l][1] = variant_IDs[(var[2],var[3]+1)]+1
    l+=1
print("hic_minor allele links: "+str(hic_minor))
l = 0
local_minor = 0
for var in local_var_counts:
    for i in range(4): dat['local_linkcounts'][l][i] = local_var_counts[var][i]
    gt1 = 0
    for i in range(4):
        if dat['local_linkcounts'][l][i] > 1: gt1+=1
    if gt1 > 1:
        local_minor += 1
    dat['local_linksites'][l][0] = variant_IDs[(var[0],var[1]+1)]+1
    dat['local_linksites'][l][1] = variant_IDs[(var[2],var[3]+1)]+1
    l+=1
print("local_minor allele links: "+str(local_minor))

dat['subsample'] = 10

with open(target_contig_id+'.standat.json', 'w') as json_file:
    json.dump(dat, json_file)

#l = 0
#for vID in vID_lookup:
#    dat['site_map'][l][0] = vID_lookup[l][0]
#    dat['site_map'][l][1] = vID_lookup[l][1]
#    l+=1

exit(0)

# find all contig ends that have two outgoing edges
segmap = {}
segseqlens = {}
rev = {'+':'-','-':'+'}
for seg in gfa.segments:
    segmap[seg.name+'+'] = set()
    segmap[seg.name+'-'] = set()
    segseqlens[seg.name]=len(seg.sequence)


for edge in gfa.edges:
#    segmap[(edge.from_segment.name,edge.from_orient)].append((edge.to_segment.name,edge.to_orient))
#    segmap[(edge.to_segment.name,edge.to_orient)].append((edge.from_segment.name,edge.from_orient))
    if segseqlens[edge.from_segment.name] < MIN_SINGLECOPY_LENGTH and segseqlens[edge.to_segment.name] < MIN_SINGLECOPY_LENGTH: continue
    segmap[edge.from_segment.name+edge.from_orient].add(edge.to_segment.name+edge.to_orient)
    segmap[edge.to_segment.name+rev[edge.to_orient]].add(edge.from_segment.name+rev[edge.from_orient])

# print(segmap)

G = nx.Graph()
for end in segmap:
    for dest_i in segmap[end]:
        for dest_j in segmap[end]:
            G.add_edge(dest_i, dest_j)
            if dest_i == ('edge_2668-') or dest_j == ('edge_2668-') or \
            dest_i == ('edge_2669-') or dest_j == ('edge_2669-') or \
            dest_i == ('edge_3789-') or dest_j == ('edge_3789-'):
                print(str(dest_i) + " linked to " + str(dest_j) + ' source end: '+ str(end))

components = nx.connected_components(G)
for comp in components:
    if len(comp)>5:
        print(comp)
        incoming = 0
        for seg in segmap:
            if len(segmap[seg] & comp) > 0:
                incoming += 1
        print("Incoming: "+str(incoming))



exit(0)

bi_seg_map = {}
bi_seg_counts = {}
for seg in gfa.segments:
    for o in ['-','+']:
        if len(segmap[(seg.name,o)]) == 2:
            if segmap[(seg.name,o)][0] in bi_seg_map or segmap[(seg.name,o)][1] in bi_seg_map:
                print("by golly this graph is complicated")
                print(segmap[(seg.name,o)][0])
                print(segmap[(seg.name,o)][1])
            bi_seg_map[segmap[(seg.name,o)][0]] = segmap[(seg.name,o)][1]
            bi_seg_map[segmap[(seg.name,o)][1]] = segmap[(seg.name,o)][0]
        if len(segmap[(seg.name,o)]) > 2:
            print("Holy mackerel: "+seg.name)
            print(segmap[(seg.name,o)])


ref_lens = {}
read_links = {}
read_link_diffs = {}
read_refs = {}
bam = pysam.AlignmentFile(sys.argv[3], "rb" )
for ref in bam.references:
    ref_lens[ref] = bam.get_reference_length(ref)

for read in bam.fetch():
    orient = ''
    if read.reference_pos < CONTIG_MARGIN:
        orient = '-'
    if read.reference_pos > ref_lens[read.reference_name] - CONTIG_MARGIN:
        orient = '+'
    if orient == '': continue
    if (read.reference_name,orient) not in bi_seg_map: continue
