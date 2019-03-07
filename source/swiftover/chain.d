module swiftover.chain;

import std.algorithm : map;
import std.algorithm.comparison : max, min;
import std.array : appender, array;
import std.algorithm : splitter;
import std.ascii : isWhite, newline;
import std.conv : to;
import std.file : exists, FileException;
import std.format;
import std.range;
import std.stdio;

import intervaltree.splaytree;

import dhtslib.htslib.hts_log;

/// ASCII valued to speed assignment and cmp
enum STRAND : ubyte
{
    PLUS = '+',
    MINUS = '-'
}

/// Zero-based, half-open
struct ChainInterval
{
    string contig;  /// whatever
    int start;      /// whatever
    int end;        /// whatever
    STRAND strand;  /// + or -
    bool invertStrand;  /// whether this interval, relative to another, has opposite-strandedness

    /// This constructor necessary to allow construction compat with SplayTree::BasicInterval
    this(int start, int end)
    {
        this.start = start;
        this.end = end;
    }
    /// This constructor, even though matches data layout/default ctor, necessary due to above
    this(string contig, int start, int end)
    {
        this.contig = contig;
        this.start = start;
        this.end = end;
    }
    /// ditto
    this(string contig, int start, int end, STRAND strand, bool invert = false)
    {
        this.contig = contig;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.invertStrand = invert;
    }


    string toString() const
    {
        return format("%s:%d-%d", this.contig, this.start, this.end);
    }

    invariant
    {
        assert(this.start <= this.end, "start not <= end");
    }
}

/// Represents a mapping from ChainInterval -> ChainInterval
/// Beacuse this is what's stored in RBTree nodes,
/// should also implement Interface Interval, but delegate to member 'query'
/// (tree sorted by query, not target)
struct ChainLink
{
    // target and query intervals in 1:1 bijective relationship
    ChainInterval target;   /// Target (reference)
    ChainInterval query;    /// Query (destination)

    bool invertStrand;      /// true iff target.strand != query.strand;

    // To work with the interval tree, and overlap functions,
    // which need access to "start" and "end"
    alias target this;
    //alias start = target.start;   // fails from outside: Error: need this for start of type int
    //alias end = target.end;       // fails from outside: Error: need this for end of type int

    int delta;  /// fixed offset from target.start->query.start  == query.start - target.start (same for end as intervals must be same len)

    /// Overload <, <=, >, >= for ChainLink/ChainLin; compare query
	@nogc int opCmp(ref const ChainLink other) const nothrow
	{
		// Intervals from different contigs are incomparable
        // TODO: or are they???
		assert(this.target.contig == other.target.contig);	// formerly return 0, but implies equality, leads to "duplicate insert" bug when using std.container.rbtree
		
		if (this.start < other.start) return -1;
		else if(this.start > other.start) return 1;
		else if(this.start == other.start && this.end < other.end) return -1;	// comes third as should be less common
		else if(this.start == other.start && this.end > other.end) return 1;
		else return 0;	// would be reached in case of equality (although we do not expect)
	}
	/// Overload <, <=, >, >= for ChainLink/int 
	@nogc int opCmp(const int x) const nothrow
	{
		if (this.start < x) return -1;
		else if (this.start > x) return 1;
		else return 0;
	}

	// TODO print strand?
    string toString() const
	{
		return format("%s:%d-%d → %s:%d-%d",
                    this.target.contig, this.target.start, this.target.end,
                    this.query.contig, this.query.start, this.query.end);
	}

    invariant
    {
        // TODO: should we really allow start == end?
        assert(this.target.start <= this.target.end);
        assert(this.query.start <= this.query.end);
        assert((this.query.end - this.query.start) == (this.target.end - this.target.start),
            "ChainLink intervals differ in length");
        assert(this.delta == (this.query.start - this.target.start));

        if (this.target.strand != this.query.strand) assert(this.invertStrand, "invertStrand error");
        else assert(!this.invertStrand, "invertStrand error");
    }
}
unittest
{
    const auto c1 = ChainLink( ChainInterval(1000, 2000), ChainInterval(10_000, 11_000), 9000);
    const auto c2 = ChainLink( ChainInterval(1000, 3000), ChainInterval(10_000, 12_000), 9000);
    const auto c3 = ChainLink( ChainInterval(1500, 2500), ChainInterval(10_500, 11_500), 9000);

    assert(c1 < c2, "ChainLink opCmp problem");
    assert(c2 < c3, "ChainLink opCmp problem");
    assert(c3 < 1501, "ChainLink opCmp problem");
    assert(c1 > 999, "ChainLink opCmp problem");
}

/// The Chain type represents a UCSC chain object, including all
/// the fields from the header line and each block of mappings
/// for that chain.
/// 
/// Specifications: https://genome.ucsc.edu/goldenpath/help/chain.html
///
/// BEWARE:
/// Note that in the above documentation, UCSC calls the starting build
/// "target"* or "reference", and the liftover-to-build "query".
///
/// This is the opposite of what I would like, given that in our application
/// we will query an interval in the reference to obtain an interval in
/// a targeted destination genome build.
///
/// *actually they don't say "target", but tName, tSize, tStrand, etc.
struct Chain
{
	long score;         /// Chain score

	string targetName;  /// source build contig
	int targetSize;     /// source build contig length
	char targetStrand;	/// +,-
	int targetStart;
	int targetEnd;
	
	string queryName;   /// destination build contig
	int querySize;      /// destination build contig length
	char queryStrand;   /// +,-
	int queryStart;
	int queryEnd;

	int id;             /// chain id

    bool invertStrand;  /// whether the strand in target and query differ

	ChainLink*[] links;  /// query and target intervals in 1:1 bijective relationship

    static auto hfields = appender!(char[][]);  /// header fields; statically allocated to save GC allocations
    static auto dfields = appender!(int[]);       /// alignment data fields; statically allocated to save GC allocations

    /// Construct Chain object from a header and range of lines comprising the links or blocks
	this(R)(R lines)
    if (isInputRange!R)
	{
        this.links.reserve(8192);

        // Example chain header line: 
		// chain 20851231461 chr1 249250621 + 10000 249240621 chr1 248956422 + 10000 248946422 2
		// assumes no errors in chain line
        this.hfields.clear;
        this.hfields.put(lines.front().dup.splitter(" "));
        assert(this.hfields.data[0] == "chain");
        this.score = this.hfields.data[1].to!long;
        
        this.targetName = this.hfields.data[2].idup;
        this.targetSize = this.hfields.data[3].to!int;
        this.targetStrand = this.hfields.data[4][0];
        this.targetStart = this.hfields.data[5].to!int;
        this.targetEnd = this.hfields.data[6].to!int;

        this.queryName = this.hfields.data[7].idup;
        this.querySize = this.hfields.data[8].to!int;
        this.queryStrand = this.hfields.data[9][0];
        this.queryStart = this.hfields.data[10].to!int;
        this.queryEnd = this.hfields.data[11].to!int;

        this.id = this.hfields.data[12].to!int;
        
        if (this.targetStrand != this.queryStrand) this.invertStrand = true;

        /*  Alignment data lines:
            size dt dq
            size -- the size of the ungapped alignment
            dt -- the difference between the end of this block and the beginning of the next block (reference sequence)
            dq -- the difference between the end of this block and the beginning of the next block (query sequence)
            NOTE: The last line of the alignment section contains only one number: the ungapped alignment size of the last block.
        */
        auto tFrom = this.targetStart;   // accum current coordinate, target
        auto qFrom = this.queryStart;    // accum current coordinate, query
        bool done;
        foreach(line; lines.dropOne())
        {
            if (line.length == 0) continue; // blank lines end a chain block

            this.dfields.clear();
            this.dfields.put(line.splitter.map!(x => x.to!int));   // TODO: benchmark splitter() 
            
            immutable int size = this.dfields.data[0];
            // note that dt and dq are not present in the final row of a chain

            // set up ChainLink from alignement data line
            ChainLink* link = new ChainLink;

            link.target.contig = this.targetName;
            link.target.start = tFrom;
            link.target.end = tFrom + size;
            link.target.strand = cast(STRAND) this.targetStrand;

            // "When the strand value is "-", position coordinates
            // are listed in terms of the reverse-complemented sequence."
            // (https://genome.ucsc.edu/goldenpath/help/chain.html)
            // TODO, this could be sped up by precomputing some values outside the foreach
            if(this.queryStrand == '+') {
                link.query.contig = this.queryName;
                link.query.start = qFrom;
                link.query.end = (qFrom + size);
                link.query.strand = cast(STRAND) this.queryStrand;
            } else {
                link.query.contig = this.queryName;
                link.query.end = this.querySize - qFrom;
                link.query.start = this.querySize - (qFrom + size);
                link.query.strand = cast(STRAND) this.queryStrand;
            }

            link.delta = qFrom - tFrom;

            link.invertStrand = this.invertStrand;

            if(this.dfields.data.length == 1)    // last block in chain
                done = true;
            else if(this.dfields.data.length == 3)
            {
                assert(done is false, "Malformed alignment data blocks");

                immutable int dt = this.dfields.data[1];
                immutable int dq = this.dfields.data[2];

                tFrom += (size + dt);
                qFrom += (size + dq);
            }
            else assert(0, "Unexpected length of alignment data line");

            // store in this.links
            this.links ~= link;
        }
	}

	string toString() const
	{
		return format("%s:%d-%d → %s:%d-%d :: %d links", 
			this.targetName, this.targetStart, this.targetEnd,
            this.queryName, this.queryStart, this.queryEnd,
            this.links.length);
	}

    invariant
    {
        assert(this.targetName != "");
        assert(this.targetSize > 0);
        assert(this.targetStrand == '+' || this.targetStrand == '-');
        assert(this.targetStart < this.targetEnd);

        assert(this.queryName != "");
        assert(this.querySize > 0);
        assert(this.queryStrand == '+' || this.queryStrand == '-');
        assert(this.queryStart < this.queryEnd);
    }
}
unittest
{
    string data = r"chain 267340 chrX 155270560 + 49204606 49241588 chrX 156040895 + 49338635 49547634 916
75      964     965
128     1047    1049
686     37      36
63      12      13
279     51      51
204     73      73
73      25      28
145     38      39
100     55      55
270     2613    2636
148     242     242
1292    11570   183551
304     3114    3099
55      43      43
51      5975    6008
104     3239    3227
40      182     182
1334    2239    2239
112";

    auto c = Chain(data.splitter(newline));  // TODO, need to change this out for cross-platform \n\r \n \r splitter

    assert(c.links.length == 19,
        "Failure parsing chain data blocks into ChainLinks");
    
    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_trace(__FUNCTION__, format("Chain: %s", c));
}


/// Representation of UCSC-format liftover chain file, which contains multiple alignment/liftover chains
struct ChainFile
{
    // TODO use build names
    string sourceBuild; /// original assembly, e.g. hg19 in an hg19->GRCh38 liftover
    alias queryBuild = sourceBuild;
    string destBuild; /// destination assembly, e.g. GRCh38 in an hg19->GRCh38 liftover

    private IntervalSplayTree!(ChainLink)*[string] chainsByContig; /// AA of contig:string -> Interval Tree

    /// Parse UCSC-format chain file into liftover trees (one tree per source contig)
    this(string fn)
    {
        if (!fn.exists)
            throw new FileException("File does not exist");

        // TODO: speed this pig up
        auto chainArray = fn.File.byLineCopy().array();

        long chainStart;
        long chainEnd;
        foreach(i, line; enumerate(chainArray))
        {
            if (line.length > 5 && line[0..5] == "chain")
            {
                chainEnd = i - 1;

                if (chainEnd > 0) // first iteration does not mark the end of a chain
                {
                    auto c = Chain(chainArray[chainStart..chainEnd]);
                    debug hts_log_trace(__FUNCTION__, format("Chain: %s", c));
                    
                    // Does this contig exist in the map?
                    auto tree = this.chainsByContig.require(c.targetName, new IntervalSplayTree!ChainLink);
                    
                    // Insert all intervals from the chain into the tree
                    foreach(link; c.links)
                        (*tree).insert(*link);

                }
                chainStart = i;
            }
        }
        // Don't forget the last chain
        auto c = Chain(chainArray[chainStart .. $]);
        debug hts_log_trace(__FUNCTION__, format("Final Chain: %s", c));

        // Does this contig exist in the map?
        auto tree = this.chainsByContig.require(c.targetName, new IntervalSplayTree!ChainLink);

        // Insert all intervals from the chain into the tree
        foreach(link; c.links)
            (*tree).insert(*link);

        debug { // debug-only to speed startup
            foreach(contig; this.chainsByContig.byKey) {
                hts_log_trace(__FUNCTION__, format("Contig: %s", contig));
            }
        }

        // show me what's in the last accessed tree:
        debug { // debug-only to speed startup
            while(tree.iteratorNext() !is null)
            hts_log_trace(__FUNCTION__, format("sorted tree entry: %s", *tree.cur));
        }
        // END ChainFile ctor
    }

    /** Lift coordinates from one build to another

        TODO: error handling (at cost of speed)

        Returns:    array (TODO, range) of ChainLink
                    Initially was ChainLink, then ChainInterval, then ChainLink again
                    because we need strandInvert *and* delta *and* source/dest
    */
    ChainLink[] lift(string contig, int start, int end)
    {
        auto i = BasicInterval(start, end);
        auto o = this.chainsByContig[contig].findOverlapsWith(i);   // returns Node*(s)

        // marked as debug because in hot code path
        debug foreach(x; o) hts_log_trace(__FUNCTION__, format("%s", *x));

        if (o.length == 0)  // no match
        {
            debug hts_log_trace(__FUNCTION__, "No match to interval");

            // -O3 will elide the stack alloc, thanks godbolt.org !
            ChainLink[] ret;
            return ret;
        }
        else if (o.length == 1) // one match
        {
            debug hts_log_trace(__FUNCTION__, format("Basic interval: %s | overlap interval: %s | delta %d",
                i, o.front().interval, o.front().interval.delta));

            // intersect makes the chain link comply with bounds of interval
            const auto isect = o.front().interval.intersect(i);
            //TODO here is where we would want to report truncated interval

            debug if(isect.start + o.front().interval.delta == 248_458_169)
                hts_log_debug(__FUNCTION__, format("%s", o.front().interval));

            return [isect];
        }
        else    // TODO optimize; return Range
        {
            auto isect = o.map!(x => x.interval.intersect(i));

            return array(isect);
        }
    }
}

/// return the intersection of a ChainLink with an interval
/// NOTE: this requires the first two elements of InteravalType
/// be (start, end) or have ctor(start, end)
/// TODO: error handling (at cost of speed)
pragma(inline, true)
ChainLink intersect(IntervalType)(ChainLink int1, IntervalType int2)
if (__traits(hasMember, IntervalType, "start") &&
    __traits(hasMember, IntervalType, "end"))
{
    // Trimming is "inwards", i.e., ===== -> -===-
    const auto trimmedStart = max(int1.start, int2.start);
    const auto trimmedEnd   = min(int1.end,   int2.end);

    const auto startDiff = trimmedStart - int1.start;
    const auto endDiff   = int1.end - trimmedEnd;

    // TODO this will need a small rewrite if I simplify ChainLink
    return ChainLink(
            ChainInterval(int1.target.contig, trimmedStart, trimmedEnd),
            ChainInterval(int1.query.contig,  int1.query.start + startDiff,
                                            int1.query.end - endDiff),
            int1.invertStrand,  // invertStrand
            int1.delta);        // delta
}
/*
pragma(inline, true)
@nogc nothrow
IntervalType intersect(IntervalType)(IntervalType int1, IntervalType int2)
if (__traits(hasMember, "IntervalType", "start") &&
    __traits(hasMember, "IntervalType", "end"))
{
    return IntervalType(
        max(int1.start, int2.start),
        min(int1.end, int2.end)        
    );
}
/// ditto
pragma(inline, true)
@nogc nothrow
BasicInterval intersect(IntervalType1, IntervalType2)(IntervalType1 int1, IntervalType2 int2)
if (!is(IntervalType1 == IntervalType2) &&
    __traits(hasMember, IntervalType1, "start") &&
    __traits(hasMember, IntervalType1, "end") &&
    __traits(hasMember, IntervalType2, "start") &&
    __traits(hasMember, IntervalType2, "end"))
{
    return BasicInterval(
        max(int1.start, int2.start),
        min(int1.end, int2.end)
    );
}
*/