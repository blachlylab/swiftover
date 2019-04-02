#include <stdio.h>

#include <htslib/vcf.h>

int main(void)
{
    vcfFile *infile = bcf_open("resources/vcftest.vcf", "r");
    vcfFile *outfile = bcf_open("/tmp/out.vcf", "w");

    bcf_hdr_t* infile_hdr = bcf_hdr_read(infile);
    bcf_hdr_t* outfile_hdr = bcf_hdr_dup(infile_hdr);

    /*     void bcf_hdr_remove(bcf_hdr_t *h, int type, const char *key); */
    bcf_hdr_remove(outfile_hdr, BCF_HL_CTG, "chr3");
    bcf_hdr_sync(outfile_hdr);

    bcf_hdr_append(outfile_hdr, "##contig=<ID=chr333,length=333333>");
    bcf_hdr_append(outfile_hdr, "##contig=<ID=chr3,length=123");

    bcf_hdr_write(outfile, outfile_hdr);

    bcf_close(infile);
    bcf_close(outfile);

    return 0;
}
