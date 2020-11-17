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

import dhtslib.htslib.hts_log;

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
void liftMAF(string chainfile, string infile, string outfile, string unmatched)
{
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

    hts_log_info(__FUNCTION__, "Reading MAF");
    foreach(line; fi.byLine())
    {
        if (line.length > 0 && line[0] == '#') continue;

        fields.clear();
        fields.put( line.splitter() );

        const auto numf = fields.data.length;

        //string contig = fields.data[4].idup;
        int start = fields.data[5].to!int;
        int end = fields.data[6].to!int;

        // array (TODO: range) of matches as ChainLink(s)
        auto trimmedLinks = cf.lift(fields.data[4], start, end);

        if (trimmedLinks.length == 0)
            fu.writef("%s\n", fields.data.join("\t"));  // Write to unmatched file
        
        else // One or more resulting output intervals
        {
            // multiple output intervals could be returned in a random order; sort
            // IMPLEMENTATION NOTE: occasionally there will be output intervals with same qStart
            // but different qEnd; order is undefined
            alias querySort = (x,y) => x.qStart < y.qStart;
            foreach(link; sort!querySort(trimmedLinks))
            {
                assert(
                    link.qcid < cf.contigNames.length, 
                    format("A query contig id (%d) is not present in the array of contig names (len {%d})", link.qcid, cf.contigNames.length));
                fields.data[4] = cf.contigNames[link.qcid].dup();

                start = link.qStart;
                end   = link.qEnd;

                orderStartEnd(start, end);              // if invert, start > end, so swap
                fields.data[5] = start.toChars.array;   // 67% time vs .text.dup;
                fields.data[6] = end.toChars.array;
                
                // BED col 6: strand
                if (numf >= 6) fields.data[5][0] = STRAND_TABLE[ fields.data[5][0] + link.invert ];

                fo.writef("%s\n", fields.data.join("\t"));
            }
        }
    }
}
