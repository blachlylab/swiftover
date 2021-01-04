![swiftover logo](swiftover_logo1x.png)

# Superfast Liftover

See our preprint here: (pending)

## Background and Motivation

Our initial goal was to write a liftover tool at least as fast as Jim Kent's seminal "liftOver".
Since then, Swiftover has turned into a vehicle for exploration of interval data structures as well as a vehicle to get us high-speed VCF liftover, not available in `liftOver`, and unfortunately slow in CrossMap, the current standard.

We hypothesized that specifically for sorted genome intervals, the implicit predictive caching of splay trees could outperform other tree structures in linear/sequential search workloads as often found in genomics.
Indeed, splay trees outperform the well-balanced AVL tree, which outperformed a slightly-less-well-balanced Red-Black tree.

We have also implemented an uncommon optimization: the probabilistic splay tree. Put simply, the splaying operation is randomly performed with some probability p_i on insert and p_l on lookup. When tuned properly, this can further increase lookup speed in sequential workloads significantly.

With the recent invention of Implicit Interval Trees (IITrees) by Heng Li [1] as well as other innovative structures by others, we've tested these and found them to be even faster than play trees in _some_ linear-scan liftover workloads.
For strict, simple liftover applications as for example hg19 -> GRCh38, IITrees may slightly outperform the Splay Tree. However, for some complex liftovers, like hg19 -> CanFam4, Splay Trees outperform the competition. (**NB:** This may also be related to `GC.addRange` calls making IITree safe for GC collected memory, more tests in progress)

Our `intervaltree` library [2] provides splay trees, AVL trees, and IITrees; swiftover's backing store may be selected at compile time. **For a balance of performance in all tested cases thus far, we recommend using the Splay Trees version** (see "compilation," below) of Swiftover.
 
## Installation

Precompiled linux binaries are available on the releases page. These binaries are statically linked against htslib; a system installation of htslib is not required. This is tested and known to work on Ubuntu and Fedora; for Alpine linux and others using musl _swiftover_ must be compiled from source.

Note that the precompiled binaries are not guaranteed to be up-to-date, and a from-source installation is generally preferred, if you are able.

See also [Compiling from source](#compiling-from-source)

## Quickstart

`./swiftover -t bed -c chainfile.chain -i input.bed -u unmatched.bed > output.bed`

If lifting over MAF or VCF, `-g <destination genome.fa>` is also required.

If lifting over MAF, `-b <destination genome build name>` is also required.

## Requirements

### Chain file

If you are working with human data, you can quickly grab hg19 to hg38 (UCSC-style contig naming: chr1, chrM) and GRCh37 to GRCh38 (Ensembl-style contig naming: 1, MT) by issuing `make chains`.
Chainfiles will be placed in `resources/`, and for the time being must be un-gzipped.

Obtain chain files from UCSC or Ensembl:

[http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/)

[ftp://ftp.ensembl.org/pub/assembly_mapping/](ftp://ftp.ensembl.org/pub/assembly_mapping/)

**Swiftover needs uncompressed chain files.** In the future we will add gzip reader.

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

### Step by Step Overview

1. Install D compiler; set PATH to match
2. Ensure recent htslib is installed; set `LIBRARY_PATH` so the linker 
3. Clone swiftover repository: `git clone https://github.com/blachlylab/swiftover.git`
4. `dub build -b=release-nobounds -c=<intervaltreetype>`
    * Where `<intervaltreetype>` in { `avltree`, `splaytree`, `iitree` }

#### Install D compiler

Download for your platform here: https://github.com/ldc-developers/ldc/releases

Uncompress and include the `bin/` directory in your $PATH.

For example:
```bash
wget https://github.com/ldc-developers/ldc/releases/download/v1.24.0/ldc2-1.24.0-linux-x86_64.tar.xz
tar xfvJ ldc2-1.24.0-linux-x86_64.tar.xz
export PATH=~/ldc2-1.24.0-linux-x86_64/bin/:$PATH
```

#### Ensure recent htslib is installed

Download for your platform here: https://github.com/samtools/htslib/releases

Decompress, build, and install. If installing to the default `/user/local/lib`, be sure to set env var `LIBRARY_PATH` so the linker can find it.

For example:
```bash
wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2
tar xfvj htslib-1.11.tar.bz2
cd htslib-1.11
# if you don't have zlib, libbz2, liblzma, curl, openssl dev files the next step will yell at you
# for Debian/Ubuntu, packages are respectively: zlib1g-dev, libbz2-dev, liblzma-dev, libcurl4-openssl-dev, libssl-dev
# (the latter two are not strictly required)
./configure
make
sudo make install
# now htslib is in /usr/local/lib
# Refresh the dynamic linker cache since you've installed a new library
sudo ldconfig
```

If you didn't already set `LIBRARY_PATH`, do it now to point to whereever you've installed htslib, or else the linking step of compilation may fail if the linker does not search `/usr/local/lib`.

```bash
export LIBRARY_PATH=/usr/local/lib:$LIBRARY_PATH
```

#### Clone the swiftover repository

`git clone https://github.com/blachlylab/swiftover.git`

#### Build

Build a release version:

`dub build -b=release-nobounds`

or a debug version by default, which is significantly slower, if you have crashes or need detailed logging for some reason:

`dub build`

or, select a specific interval tree type for testing or your specialized application:

`dub build -b=release-nobounds -c=<intervaltreetype>`
    * Where `<intervaltreetype>` in { `avltree`, `splaytree`, `iitree` }

You can read in more detail in a section below with regards to tradeoffs of interval tree types.

#### Run

`./swiftover -h` for help.

If you get the error `./swiftover: error while loading shared libraries: libhts.so.3: cannot open shared object file: No such file or directory`, then the dynamic library loader was not able to find htslib, either because you did not run `ldconfig`, or because despite running this, your linux installation's ld cache configuration was not set up to scan `/usr/local/lib`. In this case, you can either add this path to the permanent configuration, typically somewhere like `/etc/ld.so.conf`, or if you do not have root access, simply set the `LD_LIBRARY_PATH` environment variable:

```bash
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```

### Other build considerations

#### Selection of interval tree type

Swiftover uses the [intervaltree](https://github.com/blachlylab/intervaltree) library.

With dub, the configurations `avltree`, `splaytree`, and `iitree` are available.
Strictly ordering these according to execution speed is not possible, as they vary for different workloads.
See background and discussion, above. We suggest splaytree for best overall performance at this time.

#### Selection of compiler

DMD codegen can be poor compared to LDC and GDC, with execution too slow to compete with `liftover`.
Use LDC2 and `dub -b=release` for > 100% speedup. Additionally, as of dklib 0.1.2, DMD cannot inline
one of the functions in khash (`kh_hash_func`), which means compilation of swiftover with LDC2 or GDC is required for best performance.

#### linking to htslib

when using LDC2, or when using the GOLD linker (instead of GNU ld), you'll need to make sure that the linker can find libhts, which is often installed in `/usr/local/lib`. GOLD does not search there by default, nor does it examine `LD_LIBRARY_PATH`. It does, however, search `LIBRARY_PATH`, so add `export LIBRARY_PATH=/usr/local/lib` (or wherever you have installed htslib) to build scripts or run before dub build.

thanks to [http://jbohren.com/articles/2013-10-28-gold/](http://jbohren.com/articles/2013-10-28-gold/)

#### htslib version

Swiftover uses our htslib binding [dhtslib](https://github.com/blachlylab/dhtslib).

dhtslib currently works with htslib-1.10 and 1.11 (which use a new API compared to 1.9), so if compiling Swiftover from source, make sure that you have a system installation matching this.

## BUGS

2019-08-20 For VCF, INFO and FORMAT columns related to allele frequencies and calculations
may no longer be accurate in the destination genome build due to subtle mapping differences
but _especially_ if the reference allele has changed. Look for the `refchg` tag on those rows.

2020-11-18 For MAF, we don't output records that multiple matching output intervals in the destination coordinate space. A metric tracks and reports this at the end.

## References

[1] https://github.com/lh3/cgranges

[2] https://github.com/blachlylab/intervaltree
