# Swift Liftover

## Quickstart

`./swiftover -t bed -c chainfile.chain -i input.bed -u unmatched.bed > output.bed`

## Background and Motivation

Our goal is to be at least as fast as Jim Kent's seminal "liftOver" tool.

## Requirements

Chain file: Obtain from UCSC or Ensembl:

[http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/)

[ftp://ftp.ensembl.org/pub/assembly_mapping/](ftp://ftp.ensembl.org/pub/assembly_mapping/)

Swiftover needs uncompressed chain files.

## File Formats
### BED

All BED formats supported, but columns  beyond the first three (i.e., strand, thickStart, thickStop) are not yet updated.

### VCF

_WIP_