module swiftover.bed;

import std.algorithm : splitter;
import std.algorithm.sorting;
import std.array : appender, join, array;
import std.conv;
import std.file;
import std.range.primitives;
import std.stdio;

import swiftover.chain;

import dhtslib.htslib.hts_log;

/**
    Lookup table for strand (+, -) inversion

    '+' == 0x2B
    '-' == 0x2D

    By indexing into table the sum of
    invert ∈ {-1, +1} and char:strand
    we obtain the new strand char.

    33d == '!', to clue me in that something went wrong
*/
static char[256] STRAND_TABLE = [
//  0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33,'-', 33,'+', 33,'-', 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
];

/** Lift rows of infile to outfile using liftover chains in chainfile

/// BED is zero-based, half-open
/// https://genome.ucsc.edu/FAQ/FAQformat.html#format1
/// https://useast.ensembl.org/info/website/upload/bed.html
/// https://software.broadinstitute.org/software/igv/BED
///
    Rows with complete or partial match(es) will go to outfile.

    Rows with no match will go to unmatched.

    Implementation note: if an interval incompletely lifts over, for example,
    if only the first 80 of 100 nt have a lifted over representation (and no
    other chain picks up any of the last 20 nt), the truncated interval in 
    destination coordinates will be output, **but the unmatched portion will
    not go into the unmatched file for reasons of speed**. This behavior could
    be altered later if desired, but with speed penalty.
*/
void liftBED(string chainfile, string infile, string outfile, string unmatched)
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

    hts_log_info(__FUNCTION__, "Reading BED");
    foreach(line; fi.byLine())
    {
        if (line.length > 0 && line[0] == '#') continue;
        else if (line.length > 4 && (line[0..5] == "track" || line[0..5] == "brows")) continue;

        fields.clear();
        fields.put( line.splitter() );

        const auto numf = fields.data.length;

        // not sure if the .idup copy and extra storage is elided by compiler but I replaced 'contig's 2 instances with fields.data[0] directly
        //string contig = fields.data[0].idup;
        int start = fields.data[1].to!int;
        int end = fields.data[2].to!int;

        // array (TODO: range) of matches as ChainLink(s)
        auto trimmedLinks = cf.lift(fields.data[0], start, end);

        int thickStart;
        int thickEnd;
        bool hasThickInterval;
        if (numf >= 8) {
            thickStart = fields.data[6].to!int;
            thickEnd = fields.data[7].to!int;

            if (thickStart != thickEnd) {
                assert(thickEnd > thickStart);
                hasThickInterval = true;
                // perform single coordinate liftovers of thickStart/thickEnd
                if (!cf.liftCoordOnly(fields.data[0], thickStart) ||
                    !cf.liftCoordOnly(fields.data[0], thickEnd)) hasThickInterval = false;
                // TODO: need to somehow label these rows [those losing their thickInterval]; it would be too expensive to search for nearest thickStart/End
                // I believe this would correspond to liftOver's -fudgeThick option
                else orderStartEnd(thickStart, thickEnd);    // can't check link.invert and conditionally swap them because link.invert not available until below
            }
        }

        if (trimmedLinks.length == 0)
            fu.writef("%s\n", fields.data.join("\t"));  // Write to unmatched file
        else // One or more resulting output intervals
        {
            // multiple output intervals could be returned in a random order; sort
            // IMPLEMENTATION NOTE: occasionally there will be output intervals with same qStart but different qEnd; order is undefined
            alias querySort = (x,y) => x.qStart < y.qStart;
            foreach(link; sort!querySort(trimmedLinks))
            {
                // Transition to ChainLinks with embedded ids (not strings)
                //fields.data[0] = link.qContig.dup;
                assert(link.qcid < cf.contigNames.length, "A query contig id is not present in the array of contig names");
                fields.data[0] = cf.contigNames[link.qcid].dup();

                start = link.qStart;
                end   = link.qEnd;

                orderStartEnd(start, end);              // if invert, start > end, so swap
                fields.data[1] = start.toChars.array;   // 67% time vs .text.dup;
                fields.data[2] = end.toChars.array;
                
                // BED col 6: strand
                if (numf >= 6) fields.data[5][0] = STRAND_TABLE[ fields.data[5][0] + link.invert ];

                if (numf >= 8) {
                    // UCSC: "When there is no thick part, thickStart and thickEnd are usually set to the chromStart position."
                    if (!hasThickInterval || thickStart >= end || thickEnd <= start) {
                        fields.data[6] = fields.data[1].dup;    // (== start == link.qStart)
                        fields.data[7] = fields.data[1].dup;    // (== start == link.qStart)
                    }
                    else {
                        // There is a thick interval, either overlapping at front, fully contained, or overlapping at back
                        assert(thickStart < end);
                        assert(thickEnd > start);
                        
                        if (thickStart < start) {
                            fields.data[6] = start.toChars.array;
                        } else {
                            fields.data[6] = thickStart.toChars.array;
                        }
                        if (thickEnd > end) {
                            fields.data[7] = end.toChars.array;
                        } else {
                            fields.data[7] = thickEnd.toChars.array;
                        }
                    }

                } // end if (numf >= 8)
                fo.writef("%s\n", fields.data.join("\t"));
            }
        }
    }
}
