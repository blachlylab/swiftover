module swiftover.maf;

import std.algorithm : splitter;
import std.algorithm.sorting;
import std.array : appender, join, array;
import std.conv;
import std.file;
import std.format;
import std.range : enumerate;
import std.range.primitives;
import std.stdio;

import swiftover.chain;

import htslib.hts_log;
import dhtslib.faidx;

/// Nucleotide complementarity table
/// (Not including IUPAC ambiguity codes)
/// '-' (specific to MAF) left as '-'
static immutable(char)[128] NT_COMP_TABLE = [
//  0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  // 0
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  // 1
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,'-',  0,  0,  // 2
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  // 3
    0,'T',  0,'G',  0,  0,  0,'C',  0,  0,  0,  0,  0,  0,  0,  0,  // 4
    0,  0,  0,  0,'A',  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  // 5
    0,'t',  0,'g',  0,  0,  0,'c',  0,  0,  0,  0,  0,  0,  0,  0,  // 6
    0,  0,  0,  0,'a',  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]; // 7

/// Reverse complement a sequence; needed for liftover strand flips
string reverse_complement(const(char)[] allele)
{
    assert(allele.length > 0);

    // Special case "NA" which appears in MAF files
    if (allele == "NA")
        return "NA";

    // Special case len 1
    if (allele.length == 1) {
        return [NT_COMP_TABLE[ allele[0] ]];
    }
    else {
        char[] rcomp;
        rcomp.length = allele.length;
        auto len = allele.length;
        foreach(i, base; allele.enumerate(1)) {
            rcomp[len - i] = NT_COMP_TABLE[ base ];
        }
        return rcomp.idup;
    }
}

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
                Option 1: Lift range (use both), warn if length2 != length1
                Option 2: Lift one, do arithmetic for the end coord

        Strand:
            Included as col. 8, but specifications @ GDC page indicate that
            strand is always reported as '+'
                This also reminds me that we need to be careful to lift
                the end coordinate since it could reverse orientation
                (i.e., must choose Option 1 above)

        May contain header "#version 2.x.y" (currently 2.4.1?)
        Second line is apparently ALSO a header, but doesn't contain # sigil
        It is a column list (Hugo_Symbol\tEntrez_Gene_Id\tCenter\t...)

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
    hts_log_warning(__FUNCTION__, "Only fields build, chr, start, end, ref/tumor1/tumor2 are updated.");
    hts_log_warning(__FUNCTION__, "A reference allele change might invalidate other fields.");

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

    int nmatched, nunmatched, nmultiple, nrefchg, nneither;

    hts_log_info(__FUNCTION__, "Reading MAF");
    // Consume header lines; write to output and unmatched
    string version_header = fi.readln();
    string fields_header = fi.readln();
    fo.write(version_header);
    fu.write(version_header);
    fo.write(fields_header);
    fu.write(fields_header);
    foreach(line; fi.byLine())
    {
        fields.clear();
        fields.put( line.splitter() );

        const auto numf = fields.data.length;

        long start = fields.data[MAF.Start].to!long;
        long end = fields.data[MAF.End].to!long;

        // MAF is 1-based, closed (i.e., 1 base is specified s.t. start==end)
        // Chain files are zbho
        assert(start <= end);
        start = start - 1;

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

            // (code below referencing `link` copied from BED where we use foreach(link; trimmedLinks)
            auto link = trimmedLinks[0];
            assert(link.qcid < cf.contigNames.length, 
                    format("A query contig id (%d) is not present in the array of contig names " ~
                            "(len {%d})", link.qcid, cf.contigNames.length));
            fields.data[MAF.Contig] = cf.contigNames[link.qcid].dup();

            start = link.qStart;
            end   = link.qEnd;

            orderStartEnd(start, end);              // if invert, start > end, so swap
            start++;                                // Revert to one based, closed coords
            fields.data[MAF.Start] = start.toChars.array;   // 67% time vs .text.dup;
            fields.data[MAF.End] = end.toChars.array;

            // If there was a strand flip, update tumor alleles
            // (Previously, this was done only in case of a ref change, but in theory there could be
            // a strand flip and NO reference change, not negating the need to rcomp other alleles)
            // Note that because of the checks in the REF change block below, we have to do this first.
            if (link.invert < 0)
            {
                fields.data[MAF.Tumor_Seq_Allele1] = reverse_complement(fields.data[MAF.Tumor_Seq_Allele1]).dup;
                fields.data[MAF.Tumor_Seq_Allele2] = reverse_complement(fields.data[MAF.Tumor_Seq_Allele2]).dup;
            }

            /* Check Reference Allele

                Obviously, this needs to be done after contig/start/end have
                been lifted into destination coordinates.

                In addition, MAF requires some special handling. In particular,
                the REF field of records of type insertion is listed as `-`, rather
                than the single nucleotide at the position, as in VCF. For type
                deletion, it includes the deleted sequence.

                Unfortunately in the case of an insertion, without an actual REF
                allele (from source build) we cannot know if there has been a
                reference allele change. The GDC MAF format does include a "context"
                field, but this may not be present in all MAFs. Per the spec,
                "Novel inserted sequence for insertion does not include flanking
                reference bases."

                In the case of deletion, the ref allele is given and `-` symbol
                represents the variant. AFAICT, there is some inconsistency, at least
                in the MSKCC cancerhotspots v2 MAF, whereby sometimes deletions'
                tumor alleles are shown as (NA, -), but sometimes (-, NA), and
                sometimes (REF, -) but othertimes not. Unclear.

                Plan:
                SNP/SNV: Check reference allele, track if changed; warn if the new
                    ref allele matches NEITHER of the tumor detected alleles, which
                    is certainly possible but may be notable (TODO: record or report
                    these locations)

                Insertion: Skip reference allele check
                    NOTE: It is possible that a tumor "somatic" insertion is in fact
                    the new reference; this could in theory be detected if the liftover
                    coordinates were of different length (i.e. > 1) than the source coords.
                    We can special case this later if necessary

                Deletion: TODO until we can figure out whether our reference MAF
                    is actually systematic (see above) or not, and what is the pattern
                    with respect to how the tumor seq allele1/2 are represented
            */
            auto newRefAllele = (fields.data[MAF.Ref_Allele] == "-") ? "-" :
                fa.fetchSequence!(CoordSystem.obc)(fields.data[MAF.Contig].idup, start, end);
            if (fields.data[MAF.Ref_Allele] != newRefAllele)
            {
                nrefchg++;

                // flag if the new reference matches NEITHER OF the somatic calls
                // (unless it's a DEL which would show (NA, -)
                // This typically occurs when tumor somatic calls match (homozygous/LOH)
                // NOTE: Tumor alleles should have already been revcomp'd in case of strand flip
                if (fields.data[MAF.Variant_Type] != "DEL" &&
                    newRefAllele != fields.data[MAF.Tumor_Seq_Allele1] && 
                    newRefAllele != fields.data[MAF.Tumor_Seq_Allele2]) {

                    nneither++;
                    debug hts_log_warning(__FUNCTION__,
                        format("New REF matches NEITHER tumor allele: %s %d-%d %s -> %s (%s | %s)",
                            fields.data[MAF.HUGO_Symbol],
                            link.tStart,
                            link.tEnd,
                            fields.data[MAF.Ref_Allele],
                            newRefAllele,
                            fields.data[MAF.Tumor_Seq_Allele1],
                            fields.data[MAF.Tumor_Seq_Allele2]));
                }

                // update REF allele
                fields.data[MAF.Ref_Allele] = newRefAllele.dup;
            }


            // Finally, update the genome build name
            fields.data[MAF.NCBI_Build] = genomebuild.dup;
            fo.writef("%s\n", fields.data.join("\t"));
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

    hts_log_info(__FUNCTION__, format("Matched %d records (%d multiply, skipped) and %d unmatched",
                                nmatched, nmultiple, nunmatched));
    hts_log_info(__FUNCTION__, format("%d REF allele changes (%d not matching either tumor allele)",
                                nrefchg, nneither));

    if (nmatched == 0)
        hts_log_warning(__FUNCTION__,
            "Could you have a chromosome name mismatch between chain file and MAF?");

}
