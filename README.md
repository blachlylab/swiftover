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

Swiftover needs uncompressed chain files.

## File Formats
### BED

All BED formats supported, including column 6 (strand) and columns 7-8 (thickStart/thickEnd).

*CAVEATS:* swiftover does not join intervals that are discontiguous
in the destination coordinates, whereas UCSC liftOver does.

### VCF

_WIP_