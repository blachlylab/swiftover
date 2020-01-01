module swiftover.main;

import std.file;
import std.getopt;
import std.path;
import std.stdio;

import swiftover.bed;
import swiftover.vcf;

import dhtslib.htslib.hts_log;

version(avl) enum treeTypeString = " (version: AVL trees)";
version(splay) enum treeTypeString = " (version: splay trees)";
version(iitree) enum treeTypeString = " (version: cgranges/IITree)";
//else enum treeTypeString = " (version: ???)";

int main(string[] args)
{
    string fileType;
    string chainfile;
    string genomefile;
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
        "g|genome", "Genome (destination build; req. for VCF)", &genomefile,
        "i|infile", "Input file; - or omit for stdin", &infile,
        "o|outfile", "Output file; - or omit for stdout", &outfile,
        "u|unmatched", "Unmatched output file", &unmatched,
    );
    }
    catch (GetOptException e)
    {
        // TODO WTF does this not work?
        //usage.helpWanted = true;
        //defaultGetoptPrinter("swift liftover", usage.options);

        // Workaround: https://forum.dlang.org/post/smqkbkthzfvbvygzfuiz@forum.dlang.org
        auto helpflag = ["swiftover executable", "-h"];
        usage = getopt(
            helpflag,
            std.getopt.config.required,
            "t|type", "File type: bed|vcf", &fileType,
            std.getopt.config.required,
            "c|chainfile", "UCSC-format chain file", &chainfile,
            "g|genome", "Genome (destination build; req. for VCF)", &genomefile,
            "i|infile", "Input file; - or omit for stdin", &infile,
            "o|outfile", "Output file; - or omit for stdout", &outfile,
            "u|unmatched", "Unmatched output file", &unmatched,
        );
        defaultGetoptPrinter("🚀 swift liftover" ~ treeTypeString,
            usage.options);

        return 1;
    }

    if (usage.helpWanted)
    {
        defaultGetoptPrinter("🚀 swift liftover" ~ treeTypeString,
            usage.options);
        return 1;
    }

    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
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

    if (unmatched == "")
        throw new Exception("<unmatched> output file required");
    
    hts_log_info(__FUNCTION__, "🚀 swift liftover" ~ treeTypeString);
    switch(fileType)
    {
        case "bed":
            liftBED(chainfile, infile, outfile, unmatched);
            break;
        case "vcf":
            if (genomefile == "")
                throw new Exception("Genome FASTA (for the destination build)required.");
            liftVCF(chainfile, genomefile, infile, outfile, unmatched);
            break;
        default:
            throw new Exception("Unknown file type. Use \"bed\" or \"vcf\".");
    }

    version(instrument)
    {
        import std.format : format;
        import std.algorithm : sort;
        import std.algorithm.comparison : min, max;
        import std.algorithm.searching : minElement, maxElement;
        import std.algorithm.iteration : mean;
        
        double median(T)(T[] nums) pure nothrow {
            nums.sort();
            if (nums.length & 1)
                return nums[$ / 2];
            else
                return (nums[$ / 2 - 1] + nums[$ / 2]) / 2.0;
        }

        version(avl)
        {
            import intervaltree.avltree : _avltree_visited;
            auto n = _splaytree_visited.length;
            auto mind = minElement(_splaytree_visited);
            auto maxd = maxElement(_splaytree_visited);
            auto meand = mean(_splaytree_visited);
            auto mediand = median(_splaytree_visited);
        }
        version(splay)
        {
            import intervaltree.splaytree : _splaytree_visited;
            auto n = _splaytree_visited.length;
            auto mind = minElement(_splaytree_visited);
            auto maxd = maxElement(_splaytree_visited);
            auto meand = mean(_splaytree_visited);
            auto mediand = median(_splaytree_visited);
        }
        hts_log_info(__FUNCTION__, format("Tree statistics: N=%d (%d -- %d) mu=%f median %f", n, mind, maxd, meand, mediand) );

    }

    return 0;
}


