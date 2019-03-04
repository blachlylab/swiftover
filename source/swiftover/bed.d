module swiftover.bed;

import std.algorithm : splitter;
import std.array : appender, join;
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
    bool:invertStrand and char:strand
    we obtain the new strand char.

    33d == '!', to clue me in that something went wrong
*/
static char[256] STRAND_TABLE = [
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, '+', '-', '-', '+', 33, 
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

        // array (TODO: range) of matches as ChainInterval(s)
        auto ret = cf.lift(contig, start, end);

        if (ret.length == 0)
            fu.writef("%s\n", fields.data.join("\t"));
        else if (ret.length == 1)
        {
            fields.data[0] = ret.front().contig.dup;    
            fields.data[1] = ret.front().start.text.dup;    // TODO benchmark vs .toChars.array
            fields.data[2] = ret.front().end.text.dup;  // TODO benchmark vs .toChars.array
            
            //if (numf >= 6) fields.data[5][0] = ret.front().strand;
            ////if (numf >= 6) fields.data[5][0] = STRAND_TABLE[ ret.front().strand + ret.front().invertStrand ];

            fo.writef("%s\n", fields.data.join("\t"));
        }
        else    // unrolled to skip loop setup in case ret.len == 1
        {
            foreach(ci; ret)    //chaininterval in return
            {
                fields.data[0] = ci.contig.dup;    
                fields.data[1] = ci.start.text.dup;    // TODO benchmark vs .toChars.array
                fields.data[2] = ci.end.text.dup;  // TODO benchmark vs .toChars.array
                
                //if (numf >= 6) fields.data[5][0] = ret.front().strand;
                ////if (numf >= 6) fields.data[5][0] = STRAND_TABLE[ ret.front().strand + ret.front().invertStrand ];
                
                fo.writef("%s\n", fields.data.join("\t"));
            }
        }
    }
}
