![swiftover logo](swiftover_logo1x.png)

# Superfast Liftover

See our preprint here: (pending)

## Background and Motivation

Our initial goal was to be at least as fast as Jim Kent's seminal "liftOver" tool.
Since then, Swiftover has turned into a vehicle for exploration of interval data structures as well as a vehicle to get us VCF liftover, not available in `liftOver`, and unfortunately slow in CrossMap, the current standard.

We hypothesize that specifically for sorted genome intervals, the implicit predictive caching of splay trees outperforms other tree structures in linear/sequential search workloads as often found in genomics.
Indeed, splay trees outperform the well-balanced AVL tree, which outperformed a slightly-less-well-balanced Red-Black tree.

With the recent invention of Iplicit Interval Trees (IITrees) by Heng Li [1] and similar structures by others, we've tested these and found them to be even faster than play trees in some linear-scan liftover workloads.
For strict, simple liftover applications as for example hg19 -> GRCh38, IITrees may slightly outperform the Splay Tree. However, for some complex liftovers, like hg19 -> CanFam4, Splay Trees outperform the competition. (**NB:** This may be related to `GC.addRange` calls making IITree safe for GC collected memory, more tests in progress)

Thus for a balance of performance in all tested cases thus far, we recommend using the Splay Trees version (see compilation, below) of Swiftover.
 
Interestingly, when nodes need to be added and removed (as we might need to when constructing graph genome structures) a traditional tree structure may be preferred due to the IITree's need to be entirely reindexed after each insert/delete operation. Future studies (i.e., benchmarks) are needed.

## Installation

Precompiled linux binaries are available on the releases page. These binaries are statically linked against htslib; a system installation of htslib is not required. This is tested and known to work on Ubuntu and Fedora; for Alpine linux and others using musl _swiftover_ must be compiled from source.

See also [Compiling from source](#compiling-from-source)

## Quickstart

`./swiftover -t bed -c chainfile.chain -i input.bed -u unmatched.bed > output.bed`

If lifting over MAF or VCF, `-g <destination genome.fa>` is also required`

If lifting over MAF, `-b <destination genome build name>` is also required.

## Requirements

### Chain file

If you are working with human data, you can quickly grab hg19 to hg38 (UCSC-style contig naming: chr1, chrM) and GRCh37 to GRCh38 (Ensembl-style contig naming: 1, MT) by issuing `make chains`.
Chainfiles will be placed in `resources/`, and for the time being must be un-gzipped.

Obtain chain files from UCSC or Ensembl:

[http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/)

[ftp://ftp.ensembl.org/pub/assembly_mapping/](ftp://ftp.ensembl.org/pub/assembly_mapping/)

Swiftover needs uncompressed chain files. In the future we will add gzip reader.

### Reference genome

In MAF or VCF mode, a destination (post-liftover) reference genome FASTA file is required. If a FASTA index does not exist, one will be created the first time the genome is used.

## File Formats

It is critical that the contigs appearing in the _source_ file have an entry in the chain;
otherwise the program may terminate with `range violation`. Adding error checking/handling
for this is possible, but as the check would be run once for every row of input, it could
unnecessarily slow the liftover.

Likewise, in VCF mode, contig naming must be consistent across input VCF, chain file, and destination genome.

### BED

All BED formats supported, including column 6 (strand) and columns 7-8 (thickStart/thickEnd).

*CAVEATS:* 
swiftover does not join intervals that are discontiguous
in the destination coordinates, whereas UCSC liftOver does by default.
We feel that discontiguous intervals better represent to the user the relationship between source and destination sequence.

### VCF

VCF liftover works as you would expect. 😁

Keep in mind that contig names must be consistent among the chainfile, the genome, and the VCF.

Lifting a VCF file to a new genome build additionally requres a FASTA file of the new/destination genome. If it is not already faidx indexed, an index will be created automatically. If no .fai index already exists and swiftover does not have write permission to the directory containing the genome, execution will fail.

**Reference allele change:** Occasionally, the reference allele may differ even at equivalent coordinates in different genome builds. When swiftover detects this, it will update the `REF` column of the VCF record and add the tag **refchg** to the `INFO` column. These records can then be filtered by downstream tools if necessary (e.g., `bcftools view -i 'INFO/refchg=1'`)

**CAVEATS:** `INFO` and `FORMAT` column tags related to allele frequencies and calculations may no longer be accurate in the destination geneome build (due to subtle mapping differences), but _especially_ if the reference allele has changed. We will likely add cmdline flag to strip all INFO/FORMAT tags, followed later by a plugin to recalculate select values (e.g. when refchg), or scripting (e.g. Lua) capability.

### MAF

MAF support is experimental. 

MAF liftover works similar to VCF, but additionally requires a destination genome build name,
as this is a standard MAF column: `https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/`

Reference allele changes should function identically to VCF.

**BUGS**: MAF includes boht start and end coordinates. If multiple output intervals result from
a single input row (as might be seen with SVs or longer INDELs), these records are skipped
(i.e. also not written to the unmatched file). However, a metric will track and report this.
In the future we expect to deal with these, although they may require special handling with respect
to breaking up the "reference" allele.

### GFF3

_TBD_

## Compiling from source

### Install D compiler

Download for your platform here: https://github.com/ldc-developers/ldc/releases/tag/v1.23.0

Uncompress and include the `bin/` directory in your $PATH.

### Step by Step

1. Install D compiler
2. Clone swiftover repository: `git clone https://github.com/blachlylab/swiftover.git`
3. `dub build -b=release-nobounds -c=<intervaltreetype>`
    * Where `<intervaltreetype>` in { `avltree`, `splaytree`, `iitree` }

### Selection of interval tree type

Swiftover uses the [intervaltree](https://github.com/blachlylab/intervaltree) library.

With dub, the configurations `avltree`, `splaytree`, and `iitree` are available.
Strictly ordering these according to execution speed is not possible, as they vary for different workloads.
See background and discussion, above. We suggest splaytree for best overall performance at this time.

### Selection of compiler

DMD codegen can be poor compared to LDC and GDC, with execution too slow to compete with `liftover`.
Use LDC2 and `dub -b=release` for > 100% speedup. Additionally, as of dklib 0.1.2, DMD cannot inline
one of the functions in khash (`kh_hash_func`), which means compilation of swiftover with LDC2 or GDC is required for best performance.

### linking to htslib

when using LDC2, or when using the GOLD linker (instead of GNU ld), you'll need to make sure that the linker can find libhts, which is often installed in `/usr/local/lib`. GOLD does not search there by default, nor does it examine `LD_LIBRARY_PATH`. It does, however, search `LIBRARY_PATH`, so add `export LIBRARY_PATH=/usr/local/lib` (or wherever you have installed htslib) to build scripts or run before dub build.

thanks to [http://jbohren.com/articles/2013-10-28-gold/](http://jbohren.com/articles/2013-10-28-gold/)

### htslib version

Swiftover uses our htslib binding [dhtslib](https://github.com/blachlylab/dhtslib).

dhtslib currently works only with htslib-1.9, so if compiling Swiftover from source, make sure that you have a system installation of htslib-1.9. When dhtslib finalizes htslib-1.10 (breaking ABI changes and new API) support, we will also update Swiftover to use the new dhtslib.

## BUGS

2019-08-20 For VCF, INFO and FORMAT columns related to allele frequencies and calculations
may no longer be accurate in the destination genome build due to subtle mapping differences
but _especially_ if the reference allele has changed. Look for the `refchg` tag on those rows.

2020-11-18 For MAF, we don't output records that multiple matching output intervals in the destination coordinate space. A metric tracks and reports this at the end.

## References

[1] https://github.com/lh3/cgranges

[2] https://github.com/blachlylab/intervaltree
