module swiftover.main;

import std.file;
import std.getopt;
import std.path;
import std.stdio;

import swiftover.bed;
import swiftover.vcf;

import dhtslib.htslib.hts_log;

int main(string[] args)
{
    string chainfile;
    string fileType;
    string infile;
    string outfile;
    string unmatched;

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
        "u|unmatched", "Unmatched output file", &unmatched,
    );
    }
    catch (GetOptException e)
    {
        // TODO WTF does this not work?
        usage.helpWanted = true;
        defaultGetoptPrinter("swift liftover", usage.options);
        return 1;
    }

    if (usage.helpWanted)
    {
        version(avl) enum treeTypeString = " (version: AVL trees)";
        else enum treeTypeString = " (version: splay trees)";

        defaultGetoptPrinter("🚀 swift liftover" ~ treeTypeString,
            usage.options);
        return 1;
    }

    debug hts_set_log_level(htsLogLevel.HTS_LOG_DEBUG);
    debug(trace) hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);

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
    if (unmatched)
    {
        auto fo = File(unmatched, "w");
        fo.close();
    }

    switch(fileType)
    {
        case "bed":
            liftBED(chainfile, infile, outfile, unmatched);
            break;
        case "vcf":
            liftVCF(chainfile, infile, outfile);
            break;
        default:
            throw new Exception("Unknown file type. Use \"bed\" or \"vcf\".");
    }

    return 0;
}


