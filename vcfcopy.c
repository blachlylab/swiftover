#include <stdio.h>

#include <htslib/vcf.h>

#define MAX_UNPACK BCF_UN_STR

int main(void)
{
    vcfFile *infile = bcf_open("resources/gnomad.chrY.vcf", "r");
    vcfFile *outfile = bcf_open("/tmp/out.vcf", "w");

    bcf_hdr_t* infile_hdr = bcf_hdr_read(infile);
    bcf_hdr_t* outfile_hdr = bcf_hdr_dup(infile_hdr);

    // uncomment below to increase frequency of segfault; likewise you can add a second header line to increase probabilty further
    bcf_hdr_append(outfile_hdr, "##contig=<ID=chrZ,length=57227415,source=resources/hg19ToHg38.over.chain>");
    
    bcf_hdr_write(outfile, outfile_hdr);

    
    bcf1_t *b = bcf_init1();
    b->max_unpack = MAX_UNPACK;

    while( bcf_read(infile, infile_hdr, b) >= 0 )
    {
	bcf1_t *bnew = bcf_dup(b);
 	int ret = bcf_unpack(bnew, MAX_UNPACK);
	
	bcf_translate(outfile_hdr, infile_hdr, bnew);

	if (bnew->d.n_flt > 1) break;

	bcf_write(outfile, outfile_hdr, bnew);

	bcf_destroy(bnew);
	bcf_empty(b);
    }
    
    bcf_close(infile);
    bcf_close(outfile);

    return 0;
}
