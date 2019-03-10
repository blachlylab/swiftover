module swiftover.bed;

import std.algorithm : splitter;
import std.array : appender, join, array;
import std.conv;
import std.file;
import std.range.primitives;
import std.stdio;

import swiftover.chain;

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
    auto cf = ChainFile(chainfile);
    
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

    foreach(line; fi.byLine())
    {
        if (line.length > 0 && line[0] == '#') continue;
        else if (line.length > 4 && (line[0..5] == "track" || line[0..5] == "brows")) continue;

        fields.clear();
        fields.put( line.splitter() );

        const auto numf = fields.data.length;

        string contig = fields.data[0].idup;
        int start = fields.data[1].to!int;
        int end = fields.data[2].to!int;

        // array (TODO: range) of matches as ChainLink(s)
        auto trimmedLinks = cf.lift(contig, start, end);

        if (trimmedLinks.length == 0)
            fu.writef("%s\n", fields.data.join("\t"));
        else if (trimmedLinks.length == 1)
        {
            int thickStartOffset, thickSize;

            fields.data[0] = trimmedLinks.front().qContig.dup;
            start = trimmedLinks.front().qStart;
            end   = trimmedLinks.front().qEnd;

            // BED col 7 (idx 6): thickStart
            // BED col 8 (idx 7): thickEnd
            // Work on thickStart(col 7)/thickEnd(col 8) comes first because we need original start/End values
            if (numf >= 8) {
                // UCSC: "When there is no thick part, thickStart and thickEnd are usually set to the chromStart position."
                if (fields.data[6] == fields.data[7]) {
                    fields.data[6] = start.toChars.array;
                    fields.data[7] = start.toChars.array;
                }
                else {
                    // Because thickStart/thickEnd must lie within the bounds [col2,col3)
                    // we can save a costly(???) extra lookup to the IntervalTree
                    thickStartOffset = fields.data[6].to!int -
                                            fields.data[1].to!int;
                    thickSize = fields.data[7].to!int -
                                            fields.data[6].to!int;

                    auto trimmedThickStart = start + trimmedLinks.front().invert * thickStartOffset;
                    auto trimmedThickEnd   = trimmedThickStart + trimmedLinks.front().invert * thickSize;
                    orderStartEnd(trimmedThickStart, trimmedThickEnd);

                    fields.data[6] = trimmedThickStart.toChars.array;
                    fields.data[7] = trimmedThickEnd.toChars.array; 
                }
            }

            orderStartEnd(start, end);              // if invert, start > end, so swap
            fields.data[1] = start.toChars.array;   // 67% time vs .text.dup;
            fields.data[2] = end.toChars.array;
            
            // BED col 6: strand
            if (numf >= 6) fields.data[5][0] = STRAND_TABLE[ fields.data[5][0] + trimmedLinks.front().invert ];

            fo.writef("%s\n", fields.data.join("\t"));
        }
        else    // unrolled to skip loop setup in case trimmedLinks.len == 1
        {
            foreach(ci; trimmedLinks)    //chaininterval in return
            {
                int thickStartOffset, thickSize;

                fields.data[0] = ci.qContig.dup;
                start = ci.qStart;
                end   = ci.qEnd;

                // BED col 7 (idx 6): thickStart
                // BED col 8 (idx 7): thickEnd
                // Work on thickStart(col 7)/thickEnd(col 8) comes first because we need original start/End values
                if (numf >= 8) {
                    // UCSC: "When there is no thick part, thickStart and thickEnd are usually set to the chromStart position."
                    if (fields.data[6] == fields.data[7]) {
                        fields.data[6] = start.toChars.array;
                        fields.data[7] = start.toChars.array;
                    }
                    else {
                        // Because thickStart/thickEnd must lie within the bounds [col2,col3)
                        // we can save a costly(???) extra lookup to the IntervalTree
                        thickStartOffset = fields.data[6].to!int -
                                                fields.data[1].to!int;
                        thickSize = fields.data[7].to!int -
                                                fields.data[6].to!int;

                        auto trimmedThickStart = start + ci.invert * thickStartOffset;
                        auto trimmedThickEnd   = trimmedThickStart + ci.invert * thickSize;
                        orderStartEnd(trimmedThickStart, trimmedThickEnd);

                        fields.data[6] = trimmedThickStart.toChars.array;
                        fields.data[7] = trimmedThickEnd.toChars.array; 
                    }
                }

                orderStartEnd(start, end);              // if invert, start > end, so swap
                fields.data[1] = start.toChars.array;   // 67% time vs .text.dup;
                fields.data[2] = end.toChars.array;
                
                // BED col 6: strand
                if (numf >= 6) fields.data[5][0] = STRAND_TABLE[ fields.data[5][0] + ci.invert ];

                fo.writef("%s\n", fields.data.join("\t"));
            }
        }
    }
}
