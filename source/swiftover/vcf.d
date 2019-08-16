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

/** lift VCF file from src->dst (target->query)

    Params:
        chainfile   UCSC format chain file 

*/
void liftVCF(
    string chainfile, string genomefile,
    string infile, string outfile, string unmatched,
    bool removeTags = false)
{
    import std.string : toStringz;
    
    auto cf = ChainFile(chainfile);

    auto fa = IndexedFastaFile(genomefile, true);
    
    // TODO need dhtslib to support stdin/stdout

    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    debug hts_set_log_level(htsLogLevel.HTS_LOG_DEBUG);
    
    VCFReader* fi;  // ptr since decl and init inside if/else scope makes unavailable later
    if (removeTags)
        fi = new VCFReader(infile, BCF_UN_FLT);    // BCF_UN_FLT, unpacks only thru FILTER (skips INFO, FORMAT...)
    else
        fi = new VCFReader(infile, BCF_UN_ALL);    //
    
    /* ouput file
        1. remove contigs
        2. Check destination contig names/lens against genome FASTA
        3. add destination contigs to output VCF
        4. add liftover program metadata
    */
    auto fo = VCFWriter(outfile, fi.vcfhdr);

    // Add new INFO flag, "refchg" -- REF allele changed in new build
    fo.addTag!"INFO"("refchg", 0, "Flag", "REF allele changed in new genome");

    //hts_log_info(__FUNCTION__, format("Removing %d contig entries from input VCF", fo.getHeader.sequences.length));
    // There is a bug in htslib-1.9, discovered by me
    // https://github.com/samtools/htslib/issues/842
    // After removing contig, we cannot re-add one with same name.
    /*foreach(seq; fo.getHeader.sequences) {
        bcf_hdr_remove(fo.getHeader.hdr, BCF_HL_CTG, toStringz(seq));
    }*/
    // After careful consideration, I am disabling removal of "original" contigs from the output VCF
    // as not doing so makes the notion of when and how to call bcf_translate extremely complicated
    //bcf_hdr_remove(fo.getHeader.hdr, CF_HL_CTG, null);  // remove all contigs

    // Add contigs from chainfile to output VCF
    int missingInGenome;
    int newContigsAdded;
    alias keysort = (x,y) => contigSort(x.key, y.key);
    foreach(kv; sort!keysort(cf.qContigSizes.byKeyValue.array)) {
        if (!fa.hasSeq(kv.key)) {
            //debug hts_log_debug(__FUNCTION__, format("❌ %s present in chainfile but not genome.", kv.key));
            missingInGenome++;
        }
        else if (fa.seqLen(kv.key) != kv.value)
            hts_log_warning(__FUNCTION__,
                format("🤔 %s: length %d in chainfile, %d in genome. Results may be suspect.",
                kv.key, kv.value, fa.seqLen(kv.key)));
        else {
            fo.addTag!"contig"(kv.key, kv.value, "source="~chainfile);   // contig id, length
            newContigsAdded++;
        }
    }
    if (newContigsAdded)
        hts_log_info(__FUNCTION__, format("Added %d contig entries from chainfile", newContigsAdded));
    if (missingInGenome)
        hts_log_warning(__FUNCTION__, format("%d contigs present in chainfile but not destination genome.", missingInGenome));
    
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

    foreach(rec; *fi)
    {
        // TODO: setting new header after call to bcf_translate is essentially mandatory,
        // need to wrap these in convenience fn in dhtslib
        bcf_translate(fo.getHeader.hdr, fi.getHeader.hdr, rec.line);
        rec.vcfheader = fo.vcfhdr;

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

            // Check reference allele
            //const auto newRefAllele = fa[rec.chrom, rec.pos .. (rec.pos + rec.refLen)];
            const auto newRefAllele = fa.fetchSequence(rec.chrom, rec.pos, rec.pos + rec.refLen);
            auto alleles = rec.allelesAsArray;
            if (alleles[0] != newRefAllele)
            {
                //hts_log_warning(__FUNCTION__, format("REF allele mismatch: %s -> %s", alleles[0], newRefAllele));
                
                // Check ALT -- TODO, if no REF allele will be range violation
                foreach(alt; alleles[1 .. $]) {
                    // TODO, run pluggable INFO field updates
                }

                // update REF allele (impt to do this after ALT)
                alleles[0] = newRefAllele;
                rec.alleles = alleles;

                // Add flag REFCHG
                //rec.addInfo("refchg", true);  // TODO doesn't work?
                bcf_update_info(fo.getHeader.hdr, rec.line, "refchg\0".ptr, null, 1, BCF_HT_FLAG);
            }
            //hts_log_debug("foreach", "pre writeRecord");
            //stderr.writeln(*rec.line);
            fo.writeRecord(rec);
            //fo.writeRecord(fo.getHeader.hdr, rec.line);
            //hts_log_debug("foreach", "post writeRecord");
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