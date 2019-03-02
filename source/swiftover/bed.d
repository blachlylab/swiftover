module swiftover.bed;

import std.algorithm : splitter;
import std.array : appender, join;
import std.conv;
import std.file;
import std.stdio;

import swiftover.chain;

/// BED is zero-based, half-open
/// https://genome.ucsc.edu/FAQ/FAQformat.html#format1
/// 
/// https://software.broadinstitute.org/software/igv/BED
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

    auto fields = appender!(char[][]);

    foreach(line; fi.byLine())
    {
        if (line.length > 4 && (line[0..5] == "track" || line[0..5] == "brows")) continue;

        fields.clear();
        fields.put( line.splitter() );

        string contig = fields.data[0].idup;
        int start = fields.data[1].to!int;
        int end = fields.data[2].to!int;

        cf.lift(contig, start, end);

        fields.data[0] = contig.dup;
        fields.data[1] = start.text.dup;        // TODO benchmark vs .toChars.array
        fields.data[2] = end.text.dup;          // TODO benchmark vs .toChars.array

        fo.writef("%s\n", fields.data.join("\t"));

    }
}
