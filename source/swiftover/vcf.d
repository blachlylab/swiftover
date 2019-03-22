module swiftover.vcf;

import std.format;

import swiftover.chain;

import dhtslib.vcf;
import dhtslib.htslib.hts_log;

///
void liftVCF(string chainfile, string genomefile, string infile, string outfile, string unmatched)
{
    import std.string : toStringz;
    import dhtslib.htslib.vcf;
    
    auto cf = ChainFile(chainfile);
    
    // TODO need dhtslib to support stdin/stdout

    hts_set_log_level(htsLogLevel.HTS_LOG_DEBUG);
    hts_log_info("liftvcf", infile);

    //auto fp = vcf_open("chry.vcf\0".ptr, "r"c.ptr);
    //hts_log_info(__FUNCTION__, "survived vcf_open");

    auto fi = VCFReader(infile);
    auto fo = VCFWriter(outfile, fi.vcfhdr);
    fo.writeHeader();
    auto fu = VCFWriter(unmatched, fi.vcfhdr);
    fu.writeHeader();

    // GC takes care of closing; TODO scope exit flush?

    
    int nmatched, nunmatched;

    foreach(rec; fi)
    {
        string contig = rec.chrom;
        auto coord = rec.pos;

        const auto nresult = cf.liftDirectly(contig, coord);

        rec.chrom = contig;
        rec.pos = coord;

        // If not liftable, write to unmatched
        if (!nresult) {
            nunmatched++;            
            fu.writeRecord(rec);
        }
        else {
            nmatched++;
            fo.writeRecord(rec);
        }
    }
    
    hts_log_info(__FUNCTION__, format("Matched %d records (%d unmatched)", nmatched, nunmatched));
    
}
