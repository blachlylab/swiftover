![swiftover logo](swiftover_logo1x.png)

# Superfast Liftover with Splay Trees

See our preprint here: 

Please star our repo and cite our manuscript! It makes a big difference in grant applications.

## Background and Motivation

Our goal is to be at least as fast as Jim Kent's seminal "liftOver" tool.
We hypothesize that specifically for sorted genome intervals,
the implicit predictive caching of splay trees outperforms other
data structures.

## Installation

Precompiled binaries are available on the releases page. These binaries are statically linked against htslib; a system installation of htslib is not required.

## Quickstart

`./swiftover -t bed -c chainfile.chain -i input.bed -u unmatched.bed > output.bed`

## Requirements

Chain file: Obtain from UCSC or Ensembl:

[http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/)

[ftp://ftp.ensembl.org/pub/assembly_mapping/](ftp://ftp.ensembl.org/pub/assembly_mapping/)

Swiftover needs uncompressed chain files. `gunzip` any files ending in `.gz`.

## File Formats
### BED

All BED formats supported, including column 6 (strand) and columns 7-8 (thickStart/thickEnd).

*CAVEATS:* swiftover does not join intervals that are discontiguous
in the destination coordinates, whereas UCSC liftOver does by default. We feel that discontiguous intervals better represent to the user the relationship between source and destination sequence.

### VCF

Lifting a VCF file to a new genome build additionally requres a FASTA file of the new/destination genome. If it is not already faidx indexed, an index will be created automatically. If no .fai index already exists and swiftover does not have write permission to the directory containing the genome, execution will fail.

**Reference allele change:** Occasionally, the reference allele may differ even at equivalent coordinates in different genome builds. When swiftover detects this, it will update the REF column of the VCF record and add the tag **refchg** to the INFO column. These records can then be filtered by downstream tools if necessary (e.g., `bcftools view -i 'INFO/refchg=1'`)

*CAVEATS:* INFO/FORMAT tags with numeric values (for example, allele counts or minor allele fractions) are not recalculated. This is most especially relevant when the reference allele changes.

## Compiling from source

DMD codegen can be poor compared to LDC and GDC, with execution too slow to compete with `liftover`.
Use LDC2 and `dub -b=release` for > 100% speedup.

**htslib:** when using LDC2, or when using the GOLD linker (instead of traditional GNU ld), you'll need to make sure
that the linker can find libhts, which is often installed in `/usr/local/lib`. GOLD does not search there
by default, nor does it examine `LD_LIBRARY_PATH`. It does, however, search `LIBRARY_PATH`, so add
`export LIBRARY_PATH=/usr/local/lib` to build scripts or run before dub build.

thanks to [http://jbohren.com/articles/2013-10-28-gold/](http://jbohren.com/articles/2013-10-28-gold/)
