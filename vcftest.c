#include <stdio.h>

#include "htslib/vcf.h"

int main(void)
{
    vcfFile *infile = bcf_open("resources/vcftest.vcf", "r");
    vcfFile *outfile = bcf_open("/tmp/out.vcf", "w");

    bcf_hdr_t* infile_hdr = bcf_hdr_read(infile);
    bcf_hdr_t* outfile_hdr = bcf_hdr_dup(infile_hdr);

    bcf_hdr_remove(outfile_hdr, BCF_HL_CTG, "chr3");
//    bcf_hdr_append(outfile_hdr, "##contig=<ID=chr333,length=333>");
    bcf_hdr_append(outfile_hdr, "##contig=<ID=chr3,length=123");

    int ret = bcf_hdr_write(outfile, outfile_hdr);

    bcf_hdr_destroy(infile_hdr);
    bcf_hdr_destroy(outfile_hdr);
    
    bcf_close(infile);
    bcf_close(outfile);

    return ret;
}
