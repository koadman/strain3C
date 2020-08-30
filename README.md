# strain3C

strain3C is a workflow to resolve strain genomes from metagenomic data. strain3C requires as input:

1. A polished long read assembly in GFA format
2. A name-sorted bamfile containing Hi-C read pairs mapped against the assembly contigs
3. A tabix-indexed set of variant calls in VCF format, made from the Hi-C reads mapped against the long read assembly

## Usage

```
./strain3C --gfa polished.gfa --bam hic.namesorted.bam --vcf hic.vcf.gz -resume
```
