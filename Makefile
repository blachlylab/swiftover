UCSCURL = http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
ENSURL  = ftp://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz

CHAINS = resources/hg19ToHg38.over.chain.gz resources/GRCh37_to_GRCh38.chain.gz
GNOMADEXOMES=gnomad.exomes.r2.1.1.sites.vcf.bgz

# TODO make urls platform specific (linux, macosx etc.)

.PHONY: chains
chains: ${CHAINS}

resources/hg19ToHg38.over.chain.gz:
	cd resources; wget ${UCSCURL}

resources/GRCh37_to_GRCh38.chain.gz:
	cd resources; wget ${ENSURL}


.PHONY: gnomad
gnomad:	resources/${GNOMADEXOMES}

resources/${GNOMADEXOMES}:
	cd resources; wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz

utils: liftOver liftOverMerge liftUp
	chmod +x liftOver liftOverMerge liftUp

liftOver:
	wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver

liftOverMerge:
	wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOverMerge

liftUp:
	wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftUp

resources: grch38

grch38:
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.p12.genome.fa.gz
