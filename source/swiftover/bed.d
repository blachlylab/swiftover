module swiftover.bed;

/// BED is zero-based, half-open
/// https://genome.ucsc.edu/FAQ/FAQformat.html#format1
/// 
/// https://software.broadinstitute.org/software/igv/BED
void liftBED(string chainfile, string infile, string outfile)
{
    import std.algorithm : splitter;
    import std.array : appender, join;

    File fi;
    File fo;

    if (infile == "-" || infile == "")
        fi = stdin;
    else
        fi = File(infile, "r");

    if (outfile == "-" || outfile == "")
        fo = stdout;
    else
        fo = File(outfile, "w");

    auto fields = appender!(char[][]);

    foreach(line; fi.byLine())
    {
        if (line.length > 4 && (line[0..5] == "track" || line[0..5] == "brows")) continue;

        fields.clear();
        fields.put( line.splitter() );

        fields.data[0] = "chr99".dup;
        fields.data[1] = "999".dup;
        fields.data[2] = "1001".dup;

        fo.writef("%s\n", fields.data.join("\t"));

    }
}
