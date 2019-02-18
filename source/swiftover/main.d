module swiftover.main;

import std.file;
import std.getopt;
import std.path;
import std.stdio;

int main(string[] args)
{
    string chainfile;
    string fileType;
    string infile;
    string outfile;

    GetoptResult usage;

    try {
    usage = getopt(
        args,
        std.getopt.config.required,
        "t|type", "File type: bed|vcf", &fileType,
        std.getopt.config.required,
        "c|chainfile", "UCSC-format chain file", &chainfile,
        "i|infile", "Input file; - or omit for stdin", &infile,
        "o|outfile", "Output file; - or omit for stdout", &outfile,
    );
    }
    catch (GetOptException e)
    {
        defaultGetoptPrinter("swift liftover", usage.options);
        return 1;
    }

    if (usage.helpWanted)
    {
        defaultGetoptPrinter("swift liftover",
            usage.options);
        return 1;
    }

    if (!chainfile.exists)
        throw new FileException("Chainfile does not exist");
    if (infile != "-" && infile != "" && !infile.exists)
        throw new FileException("Input file does not exist");
    // test write 
    if (outfile != "-" && outfile != "")
    {
        auto fo = File(outfile, "w");
        fo.close();
    }

    switch(fileType)
    {
        case "bed":
            liftBED(chainfile, infile, outfile);
            break;
        case "vcf":
            liftVCF(chainfile, infile, outfile);
            break;
        default:
            throw new Exception("Unknown file type. Use \"bed\" or \"vcf\".");
    }

    return 0;
}

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

///
void liftVCF(string chainfile, string infile, string outfile)
{

}