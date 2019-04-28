# Swift Liftover

## Quickstart

`./swiftover -t bed -c chainfile.chain -i input.bed -u unmatched.bed > output.bed`

## Background and Motivation

Our goal is to be at least as fast as Jim Kent's seminal "liftOver" tool.
We further hypothesize that specifically for sorted genome intervals,
the implicit predictive caching of splay trees will outperform other
tree structures.

## Requirements

Chain file: Obtain from UCSC or Ensembl:

[http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/)

[ftp://ftp.ensembl.org/pub/assembly_mapping/](ftp://ftp.ensembl.org/pub/assembly_mapping/)

If you are working with human data, you can quickly grab hg19 to hg38
(UCSC-style contig naming: chr1, chrM) and GRCh37 to GRCh38
(Ensembl-style contig naming: 1, MT) by issuing `make chains`.
Chainfiles will be placed in `resources/`, and for the time being need to be un-gzipped.

Swiftover needs uncompressed chain files. TODO: Will add gzip reader.

## File Formats

It is critical that the contigs appearing in the _source_ file have an entry in the chain;
otherwise the program will terminate with `range violation`. Adding error checking/handling
for this is possible, but as the check would be run once for every row of input, it could
unnecessarily slow the liftover.

Likewise, in VCF mode, contig naming must be consistent across input VCF, chain file, and
destination genome.

### BED

All BED formats supported, including column 6 (strand) and columns 7-8 (thickStart/thickEnd).

*CAVEATS:* swiftover does not join intervals that are discontiguous
in the destination coordinates, whereas UCSC liftOver does.

### VCF

VCF liftover works as you would expect. 😁

An extra INFO column tag `refchg` is added when the reference allele changes between the
source and destination genomes.

*CAVEATS:* INFO and FORMAT column tags related to allele frequencies and calculations may
no longer be accurate in the destination geneome build (due to subtle mapping differences),
but _especially_ if the reference allele has changed. We will likely add cmdline flag to strip
all INFO/FORMAT tags, followed later by a plugin to recalculate select values (e.g. when refchg).

## Compiling from source

DMD codegen is poor, and execution is too slow. Use LDC2 and `dub -b=release` for > 100% speedup.

when using LDC2, or when using the GOLD linker (instead of traditional GNU ld), you'll need to make sure
that the linker can find libhts, which is often installed in `/usr/local/lib`. GOLD does not search there
by default, nor does it examine `LD_LIBRARY_PATH`. It does, however, search `LIBRARY_PATH`, so add
`export LIBRARY_PATH=/usr/local/lib` to build scripts or run before dub build.

thanks to [http://jbohren.com/articles/2013-10-28-gold/](http://jbohren.com/articles/2013-10-28-gold/)
