#include <stdio.h>

#include <htslib/vcf.h>

int main(void)
{
    vcfFile *infile = bcf_open("resources/vcftest.vcf", "r");
    vcfFile *outfile = bcf_open("/tmp/out.vcf", "w");

    bcf_hdr_t* infile_hdr = bcf_hdr_read(infile);
    bcf_hdr_t* outfile_hdr = bcf_hdr_dup(infile_hdr);

    bcf_hdr_write(outfile, outfile_hdr);

    bcf1_t *b = bcf_init1();
    b->max_unpack = BCF_UN_STR;
    int success = bcf_read(infile, infile_hdr, b);
    if (success < 0) {
        printf("bcf_read failure: %d\n", success);
        return 1;
    }

    bcf1_t *b2 = bcf_dup(b);
    bcf_unpack(b2, BCF_UN_STR);
    bcf_translate(outfile_hdr, infile_hdr, b2);

    bcf_hdr_write(outfile, outfile_hdr);
    bcf_write(outfile, outfile_hdr, b2);

    bcf_hdr_destroy(infile_hdr);
    bcf_hdr_destroy(outfile_hdr);
    bcf_close(infile);
    bcf_close(outfile);

    return 0;
}
