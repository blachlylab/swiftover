module swiftover.main;

import std.algorithm.searching : endsWith;
import std.file;
import std.getopt;
import std.path;
import std.stdio;

import swiftover.bed;
import swiftover.maf;
import swiftover.vcf;

import htslib.hts_log;

version(avl) enum treeTypeString = "AVL tree";
version(splay) enum treeTypeString = "splay tree";
version(iitree) enum treeTypeString = "cgranges/IITree";
//else enum treeTypeString = "(version: ???)";

version(instrument)
{
    import std.datetime : Duration;
    import std.datetime.stopwatch;
    __gshared Duration buildTime;
    __gshared Duration liftTime;
}

int main(string[] args)
{
    string fileType;
    string chainfile;
    string genomefile;
    string genomebuild;
    string infile;
    string outfile;
    string unmatched;

    GetoptResult usage;

    try {
    usage = getopt(
        args,
        std.getopt.config.required,
        "t|type", "File type: bed|maf|vcf", &fileType,
        std.getopt.config.required,
        "c|chainfile", "UCSC-format chain file", &chainfile,
        "g|genome", "Genome FASTA (destination build; req. for MAF/VCF)", &genomefile,
        "b|build", "Genome build name (destination; req. for MAF)", &genomebuild,
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
            "t|type", "File type: bed|maf|vcf", &fileType,
            std.getopt.config.required,
            "c|chainfile", "UCSC-format chain file", &chainfile,
            "g|genome", "Genome FASTA (destination build; req. for MAF/VCF)", &genomefile,
            "b|build", "Genome build name (destination; req. for MAF)", &genomebuild,
            "i|infile", "Input file; - or omit for stdin", &infile,
            "o|outfile", "Output file; - or omit for stdout", &outfile,
            "u|unmatched", "Unmatched output file", &unmatched,
        );
        defaultGetoptPrinter("🚀 swift liftover (version: " ~ treeTypeString ~ ")",
            usage.options);

        return 1;
    }

    if (usage.helpWanted)
    {
        defaultGetoptPrinter("🚀 swift liftover (version: " ~ treeTypeString ~ ")",
            usage.options);
        return 1;
    }

    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    debug hts_set_log_level(htsLogLevel.HTS_LOG_DEBUG);
    debug(trace) hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);

    if (!chainfile.exists)
        throw new FileException("Chainfile does not exist");
    try
    if (chainfile.endsWith(".gz", ".bz2", ".Z", ".zip"))
        throw new Error("Only uncompressed chainfiles are supported (for now)");
    catch (Error e) {
        stderr.writeln(e.msg);
        return -1;
    }
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
    
    hts_log_info(__FUNCTION__, "🚀 swift liftover (version: " ~ treeTypeString ~ ")");
    switch(fileType)
    {
        case "bed":
            liftBED(chainfile, infile, outfile, unmatched);
            break;
        case "maf":
            if (genomefile == "")
                throw new Exception("Genome FASTA (for the destination build) required.");
            if (genomebuild == "")
                throw new Exception("Genome build name (for the destination) required.");
            liftMAF(chainfile, genomefile, genomebuild, infile, outfile, unmatched);
            break;
        case "vcf":
            if (genomefile == "")
                throw new Exception("Genome FASTA (for the destination build) required.");
            liftVCF(chainfile, genomefile, infile, outfile, unmatched);
            break;
        default:
            throw new Exception("Unknown file type. Use \"bed\" or \"vcf\".");
    }

    // If run was instrumented, report final statistics
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
            auto n = _avltree_visited.length;
            auto mind = minElement(_avltree_visited);
            auto maxd = maxElement(_avltree_visited);
            auto meand = mean(_avltree_visited);
            auto mediand = median(_avltree_visited);
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
        version(iitree)
        {
            import core.stdc.stdlib : free;
            import intervaltree.iitree : _iitree_visited, _iitree_visited_size, _iitree_visited_capacity;
            auto _iit_visited = _iitree_visited[0 .. _iitree_visited_size]; // make D dynamic array with zero-copy :)
            auto n = _iit_visited.length;
            auto mind = minElement(_iit_visited);
            auto maxd = maxElement(_iit_visited);
            auto meand = mean(_iit_visited);
            auto mediand = median(_iit_visited);
            free(_iitree_visited);  // was malloced/realloced in cgranges.c
        }
        hts_log_info(__FUNCTION__, format("Tree statistics: N=%d (%d -- %d) mu=%f median %f", n, mind, maxd, meand, mediand) );

        hts_log_info(__FUNCTION__, format("Tree build time (msec): %d", buildTime.total!"msecs"));
        hts_log_info(__FUNCTION__, format("Tree lift  time (msec): %d", liftTime.total!"msecs"));

        stderr.writefln("%s,%d,%d,%d,%f,%f,%d,%d", treeTypeString, n, mind, maxd, meand, mediand, buildTime.total!"msecs", liftTime.total!"msecs");
    }

    return 0;
}


