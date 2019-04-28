# TODO make urls platform specific (linux, macosx etc.)

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
