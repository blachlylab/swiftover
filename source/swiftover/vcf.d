module swiftover.vcf;

import std.algorithm.sorting;
import std.array;
import std.conv;    // to!
import std.format;
import std.stdio;
import std.string;

import swiftover.chain;

import dhtslib.vcf;
import dhtslib.htslib.vcf;
import dhtslib.htslib.hts_log;
import dhtslib.faidx;

///
void liftVCF(string chainfile, string genomefile, string infile, string outfile, string unmatched)
{
    import std.string : toStringz;
    
    auto cf = ChainFile(chainfile);

    auto fa = IndexedFastaFile(genomefile, true);
    
    // TODO need dhtslib to support stdin/stdout

    hts_set_log_level(htsLogLevel.HTS_LOG_DEBUG);
    hts_log_info("liftvcf", infile);

    auto fi = VCFReader(infile);
    
    /* ouput file
        1. remove contigs
        2. Check destination contig names/lens against genome FASTA
        3. add destination contigs to output VCF
        4. add liftover program metadata
    */
    auto fo = VCFWriter(outfile, fi.vcfhdr);

    stderr.writeln("Contigs from fo (copied from fi)");
    stderr.writefln("total: %d", fo.getHeader.sequences.length);
    foreach(seq; fo.getHeader.sequences) {
        bcf_hdr_remove(fo.getHeader.hdr, BCF_HL_CTG, toStringz(seq));
    }

    alias keysort = (x,y) => contigSort(x.key, y.key);
    foreach(kv; sort!keysort(cf.qContigSizes.byKeyValue.array)) {
        if (!fa.hasSeq(kv.key))
            hts_log_warning(__FUNCTION__, format("❌ %s present in chainfile but not genome.", kv.key));
        else if (fa.seqLen(kv.key) != kv.value)
            hts_log_warning(__FUNCTION__,
                format("🤔 %s: length %d in chainfile, %d in genome. Results may be suspect.",
                kv.key, kv.value, fa.seqLen(kv.key)));
        else {
            auto ret = fo.addTag!"contig"(kv.key, kv.value, "source="~chainfile);   // contig id, length
            stderr.writeln("wrote " ~ kv.key ~ " with return value ", ret);
        }
    }
    fo.addFiledate();
    fo.addHeaderLineKV("liftoverProgram", "swiftover");
    fo.addHeaderLineKV("liftoverProgramURL", "https://github.com/blachlylab/swiftover");
    fo.addHeaderLineKV("liftoverChainfile", chainfile);
    fo.addHeaderLineKV("liftoverGenome", genomefile);

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

/// sort contigs according to custom, not strictly lexicographically
/// i.e., 1, 2, 3...10, 11, ... X, Y, MT, unlocalized/unplaced
bool contigSort(string x, string y)
{
    // "_" containing (e.g. chr11_KI270927v1_alt) should sort last
    if (x.indexOf('_') >= 0 && y.indexOf('_') == -1) return false;
    else if (x.indexOf('_') == -1 && y.indexOf('_') >= 0) return true;

    // Now either (a) both or (b) neither contain "_" (e.g. chr11_KI270927v1_alt)

    // If different numerically, sort by number, not lexicographically
    if (x.numeric != y.numeric) return x.numeric < y.numeric;
    // If same numerically, sort lexicographically
    else return x < y;
}

// Get numeric value
private 
pure @safe
int numeric(string z)
{
    if (z.length >= 3)
        if (z[0..3] == "chr") z = z[3 .. $];

    // special case X, Y, M/MT
    if (z[0] == 'X') return 100;
    else if (z[0] == 'Y') return 101;
    else if (z[0] == 'M') return 102;

    auto firstNonNum = indexOfNeither(z, "0123456789");
    if (firstNonNum == -1)
        return z.to!int;
    else if (firstNonNum == 0)  // e.g. chrUn_GL000195v1 -> Un_GL000195v1
        return 200;
    else
        return z[0 .. firstNonNum].to!int;
}