module swiftover.maf;

import std.algorithm : splitter;
import std.algorithm.sorting;
import std.array : appender, join, array;
import std.conv;
import std.file;
import std.format;
import std.range.primitives;
import std.stdio;

import swiftover.chain;

import htslib.hts_log;
import dhtslib.faidx;

// Columns index
enum MAF
{
    HUGO_Symbol,
    Entrez_Gene_id,
    Center,
    NCBI_Build,
    Contig,
    Start,
    End,
    Strand,
    Variant_Classification,
    Variant_Type,
    Ref_Allele,
    Tumor_Seq_Allele1,
    Tumor_Seq_Allele2
}
/** Lift rows of infile to outfile using liftover chains in chainfile

    MAF is a not-awesome format, especially in the era of standardized
    VCFs and standardized annotation fields. However, we are unfortunately
    stuck with it (for now) due to TCGA, GDC, and cBioPortal all using MAF.

    MAF (Mutation Annotation Format) is used by GDC
    https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/

    Its source documentation is apparently at NCI's TCGA wiki, however,
    they have the wiki read-protected requiring an NCI account
    (which must be an error?)
    https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification

    Briefly:
        1-based, inclusive
            This is the same as GTF/GFF2
            Length :: end - start + 1

        Includes both a start and end coordinate for a reference allele
        whereas most formats include only a start coordinate
            This forces some tradeoffs -- do we liftover both? There is
            potential for a gap in the destination coordinate space,
            such that the length could change.
                Option 1: Lift both, warn if length2 != length1
                Option 2: Lift one, do arithmetic for the end coord

        Strand:
            Included as col. 8, but specifications @ GDC page indicate that
            strand is always reported as '+'
                This also reminds me that we need to be careful to lift
                the end coordinate since it could reverse orientation
                (i.e., must choose Option 1 above)

        May contain header "#version 2.x.y" (currently 2.4.1?)
            
    Leftover notes from BED:
    Implementation note: if an interval incompletely lifts over, for example,
    if only the first 80 of 100 nt have a lifted over representation (and no
    other chain picks up any of the last 20 nt), the truncated interval in 
    destination coordinates will be output, **but the unmatched portion will
    not go into the unmatched file for reasons of speed**. This behavior could
    be altered later if desired, but with speed penalty.
*/
void liftMAF(
    string chainfile, string genomefile, string genomebuild,
    string infile, string outfile, string unmatched)
{
    auto fa = IndexedFastaFile(genomefile, true);
    fa.setCacheSize(1<<20); // including this improves runtime by 33%
    //fa.setThreads(1); // including this more than doubles runtime :-O

    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Reading chainfile");
    auto cf = ChainFile(chainfile);
    
    // Fee, Fi, Fo, Fum!
    File fi;
    File fo;
    File fu;

    if (infile == "-" || infile == "")
        fi = stdin;
    else
        fi = File(infile, "r");

    if (outfile == "-" || outfile == "")
        fo = stdout;
    else
        fo = File(outfile, "w");
    
    if (unmatched)
        fu = File(unmatched, "w");

    scope(exit)
    {
        fi.close();
        fo.close();
        fu.close();
    }

    auto fields = appender!(char[][]);

    int nmatched, nunmatched, nmultiple, nrefchg;

    hts_log_info(__FUNCTION__, "Reading MAF");
    foreach(line; fi.byLine())
    {
        if (line.length > 0 && line[0] == '#') continue;

        fields.clear();
        fields.put( line.splitter() );

        const auto numf = fields.data.length;

        string contig = fields.data[MAF.Contig].idup;
        long start = fields.data[MAF.Start].to!int;
        long end = fields.data[MAF.End].to!int;

        // array (TODO: range) of matches as ChainLink(s)
        auto trimmedLinks = cf.lift(fields.data[MAF.Contig], start, end);

        // If not liftable, write to unmatched file
        if (trimmedLinks.length == 0) {
            nunmatched++;
            fu.writef("%s\n", fields.data.join("\t"));
        }
        
        // One or more resulting output intervals
        else if (trimmedLinks.length == 1) {
            nmatched++;

            // Check reference allele
            auto newRefAllele = fa.fetchSequence!(CoordSystem.obc)(contig, start, end);
            if (fields.data[MAF.Ref_Allele] != newRefAllele)
            {
                nrefchg++;

                // flag if the new reference matches the somatic call
                if (newRefAllele == fields.data[MAF.Tumor_Seq_Allele1] ||
                    newRefAllele == fields.data[MAF.Tumor_Seq_Allele2])
                    hts_log_warning(__FUNCTION__, "New REF matches tumor allele");

                // update REF allele
                fields.data[MAF.Ref_Allele] = newRefAllele.dup; // newRefAllele is string
            }

            // (code below referencing `link` copied from BED where we use foreach(link; trimmedLinks)
            auto link = trimmedLinks[0];
            assert(link.qcid < cf.contigNames.length, 
                    format("A query contig id (%d) is not present in the array of contig names " ~
                            "(len {%d})", link.qcid, cf.contigNames.length));
            fields.data[MAF.Contig] = cf.contigNames[link.qcid].dup();

            start = link.qStart;
            end   = link.qEnd;

            orderStartEnd(start, end);              // if invert, start > end, so swap
            fields.data[MAF.Start] = start.toChars.array;   // 67% time vs .text.dup;
            fields.data[MAF.End] = end.toChars.array;

            // Finally, update the genome build name
            fields.data[MAF.NCBI_Build] = genomebuild.dup;
            fo.write("%s\n", fields.data.join("\t"));
        }
        else
        {
            // Due to complexity of handling potential broken-up output intervals
            // (possibly of higher likilihood compared to VCF given use of both start,end
            // versus liftDirectly(chr,start) in case of VCF) I elect to omit this for now
            nmatched++;
            nmultiple++;

            // TODO: should I write these to unmatched file? Hate to admix true unmatched from multi-match
        }
    }

    hts_log_info(__FUNCTION__, format("Matched %d records (%d multiply, skipped) and %d unmatched; %d REF allele changes",
                                nmatched, nmultiple, nunmatched, nrefchg));

    if (nmatched == 0)
        hts_log_warning(__FUNCTION__,
            "Could you have a chromosome name mismatch between chain file and MAF?");

}
