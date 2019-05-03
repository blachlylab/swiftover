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

import intervaltree;    // BasicInterval and overlaps
version(avl) import intevaltree.avltree;
version(splay) import intervaltree.splaytree;
version(iitree) import intervaltree.cgranges;

import dhtslib.htslib.hts_log;

import containers.hashmap;  // contig name -> size
import containers.unrolledlist;

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
    this(string contig, int start, int end, STRAND strand)
    {
        this.contig = contig;
        this.start = start;
        this.end = end;
        this.strand = strand;
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

/** Link between coordinate systems

    Originally, contained two structs of type ChainInterval:
    target, and query. Target represented source coordinates,
    while query was destination coordinates. This is logically
    easier to track, but together took 50 bytes, plus bool invertStrand
    and int delta in ChainLink itself meant sizeof(ChainLink) == 55.
    Goal is to get it <= 32 bytes so as to fit two ChainLink per cache line.

    Without resorting to tricks like lookup tables for contig name,
    the minimal dataset needed per ChainLink is:
    int:target.start, 
    int:target.end,
    int:query.start,
    ptr*2:query.contig,
    byte:invert (-1 or +1),
    which totals 29 bytes. Gives room to store 2* +,- if desired.

    Because this is what's stored in interval tree node,
    must implement Interface interval (if class) or at a minimum
    have members start and end; delegate or alias these to target (src),
    as tree is sorted by target (source coordinate system) not query. 
    */
struct ChainLink
{
    // target and query intervals in 1:1 bijective relationship.
    // len(target) == len(query)
    int tStart; /// Target [start   -- Zero based closed start
    int tEnd;   /// Target end)     -- Zero based open end
    int qStart; /// Query [start
    string qContig; /// Query contig
    // TODO , should will it help to repack this with string(ptr) aligned at front of struct on 16/32 byte boundary?
    
    //NB this is confusingly namd ; do not if(invert) ! perhaps beter called 'strand'
    byte invert;    /// i ∈ {-1, +1} where -1 => target/query on different strands; +1 => same strand

    // To work with the interval tree, and overlap functions,
    // which need access to "start" and "end"
    alias start = tStart;
    alias end = tEnd;

    /// Compute qEnd PRN
    pragma(inline, true)
    @nogc nothrow
    @property
    int qEnd() const
    {
        return this.qStart + (this.invert * (this.tEnd - this.tStart));
    }

    /// Compute size PRN
    pragma(inline, true)
    @nogc nothrow
    @property
    int size() const
    {
        return this.invert * (this.tEnd - this.tStart);
    }

    /// Overload <, <=, >, >= for ChainLink/ChainLin; compare query
	@safe @nogc int opCmp(ref const ChainLink other) const nothrow
	{
        /+
		// Intervals from different contigs are incomparable
        // TODO: or are they???
		assert(this.target.contig == other.target.contig);	// formerly return 0, but implies equality, leads to "duplicate insert" bug when using std.container.rbtree
		+/

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

    string toString() const
	{
        static const char[3] strand_table = [ '-', '!', '+' ];

        // Luckily, chain files always give source interval on +, so we don't need to store
		return format("%s:%d-%d(+) → %s:%d-%d(%s)",
                    "unk", this.tStart, this.tEnd,
                    this.qContig, this.qStart, this.qStart + (this.tEnd - this.tStart),
                        strand_table[this.invert + 1]);
	}

    invariant
    {
        // i ∈ {-1, +1} where -1 => target/query on different strands; +1 => same strand
        assert(this.invert == -1 || this.invert == 1, "invert EINVAL");

        // Target always with respect to + strand and thus start < end
        assert(this.tStart < this.tEnd);
        // This may not be true for QUERY if on - strand
        /+ cannot call qEnd() or .size() from invariant -- infinite loop
        if(this.invert == 1)
            assert(this.qStart < this.qEnd);
        else if(this.invert == -1)
            assert(this.qStart > this.qEnd);
        else
            assert(0);
        
        assert(this.size == (this.tEnd - this.tStart));
        +/
    }
}
unittest
{
    const auto c1 = ChainLink(1000, 2000, 10_000, 1);
    const auto c2 = ChainLink(1000, 3000, 10_000, 1);
    const auto c3 = ChainLink(1500, 2500, 10_500, -1);
    
    assert(c1 < c2, "ChainLink opCmp problem");
    assert(c2 < c3, "ChainLink opCmp problem");
    assert(c3 < 1501, "ChainLink opCmp problem");
    assert(c1 > 999, "ChainLink opCmp problem");

    assert(c1.qEnd == 11_000);
    assert(c2.qEnd == 12_000);
    assert(c3.qEnd == 9_500);   // minus strand

    assert(c1.size == 1000);
    assert(c3.size == 1000);
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
	int targetStart;    /// source target interval zero-based closed [start,
	int targetEnd;      /// source target interval zero-based open ,end)
	
	string queryName;   /// destination build contig
	int querySize;      /// destination build contig length
	char queryStrand;   /// +,-
	int queryStart;     /// If - strand, COORDINATES ARE OF REVERSE COMPLEMENT (in chain file)
	int queryEnd;       /// If - strand, COORDINATES ARE OF REVERSE COMPLEMENT (in chain file)

	int id;             /// chain id

    //bool invertStrand;  /// whether the strand in target and query differ
    byte invert;        /// i ∈ {-1, +1} where -1 => target/query on different strands; +1 => same strand

    UnrolledList!(ChainLink *) *links;  /// query and target intervals in 1:1 bijective relationship
                                        /// ptr because is @disable this(this)

    static auto hfields = appender!(char[][]);  /// header fields; statically allocated to save GC allocations
    static auto dfields = appender!(int[]);       /// alignment data fields; statically allocated to save GC allocations

    /// Construct Chain object from a header and range of lines comprising the links or blocks
	this(R)(R lines)
    if (isInputRange!R)
	{
        this.links = new UnrolledList!(ChainLink *);

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
        // "When the strand value is "-", position coordinates
        // are listed in terms of the reverse-complemented sequence."
        // (https://genome.ucsc.edu/goldenpath/help/chain.html)
        if (this.queryStrand == '+')
        {
            this.queryStart = this.hfields.data[10].to!int;
            this.queryEnd = this.hfields.data[11].to!int;
        }
        else
        {
            // in this case start > end
            this.queryStart = this.querySize - this.hfields.data[10].to!int;
            this.queryEnd = this.querySize - this.hfields.data[11].to!int;
        }

        this.id = this.hfields.data[12].to!int;
        
        //if (this.targetStrand != this.queryStrand) this.invertStrand = true;
        this.invert = (this.targetStrand == this.queryStrand) ? 1 : -1;

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
            this.dfields.put(line.splitter.map!(x => x.to!int)); 
            
            immutable int size = this.dfields.data[0];
            // note that dt and dq are not present in the final row of a chain

            // set up ChainLink from alignement data line
            ChainLink* link = new ChainLink;

            //link.target.contig = this.targetName;
            link.tStart = tFrom;
            link.tEnd = tFrom + size;
            //link.target.strand = cast(STRAND) this.targetStrand;

            link.qContig = this.queryName;
            link.qStart = qFrom;
            // link.qEnd a computed property

            link.invert = this.invert;
            /+            
            // "When the strand value is "-", position coordinates
            // are listed in terms of the reverse-complemented sequence."
            // (https://genome.ucsc.edu/goldenpath/help/chain.html)
            // TODO, this could be sped up by precomputing some values outside the foreach
            if(this.queryStrand == '+') {
                link.query.contig = this.queryName;
                link.query.start = qFrom;
                link.query.end = (qFrom + size);
                link.query.strand = cast(STRAND) this.queryStrand;
                link.delta = qFrom - tFrom;
            } else {
                link.query.contig = this.queryName;
                link.query.end = this.querySize - qFrom;
                link.query.start = this.querySize - (qFrom + size);
                link.query.strand = cast(STRAND) this.queryStrand;
                link.delta = link.query.start - tFrom;
            }
            link.invertStrand = this.invertStrand;
            +/

            if(this.dfields.data.length == 1)    // last block in chain
                done = true;
            else if(this.dfields.data.length == 3)
            {
                assert(done is false, "Malformed alignment data blocks");

                immutable int dt = this.dfields.data[1];
                immutable int dq = this.dfields.data[2];

                tFrom += (size + dt);
                qFrom += this.invert * (size + dq); // sub if on minus strand
            }
            else assert(0, "Unexpected length of alignment data line");

            // store in this.links
            *this.links ~= link;
        }
	}

	string toString() const
	{
		return format("chain %d :: %s:%d-%d → %s:%d-%d(%s) :: %d links",
            this.id, 
			this.targetName, this.targetStart, this.targetEnd,
            this.queryName, this.queryStart, this.queryEnd, this.queryStrand,
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
        if (this.queryStrand == '+')
            assert(this.queryStart < this.queryEnd);
        else
            assert(this.queryStart > this.queryEnd);
        
        if (this.targetStrand == this.queryStrand)
            assert(this.invert == 1);
        else
            assert(this.invert == -1);

        // nonempty
        assert(this.links.length > 0);
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

    auto c = Chain(data.splitter(newline));

    assert(c.links.length == 19,
        "Failure parsing chain data blocks into ChainLinks");
    
    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_trace(__FUNCTION__, format("Chain: %s", c));
}


/// Representation of UCSC-format liftover chain file, which contains multiple alignment/liftover chains
struct ChainFile
{
    /+
    string sourceBuild; /// original assembly, e.g. hg19 in an hg19->GRCh38 liftover
    alias queryBuild = sourceBuild;
    string destBuild; /// destination assembly, e.g. GRCh38 in an hg19->GRCh38 liftover
    +/

    /// AA of contig:string -> Interval Tree
    version(avl)
        private IntervalAVLTree!(ChainLink)*[string] chainsByContig;
    version(splay)
        private IntervalSplayTree!(ChainLink)*[string] chainsByContig;

    HashMap!(string,int) qContigSizes;  /// query (destination build) contigs,
                                        /// need for VCF

    /// Parse UCSC-format chain file into liftover trees (one tree per source contig)
    this(string fn)
    {
        if (!fn.exists)
            throw new FileException("File does not exist");

        // Cannot undergo static init
        this.qContigSizes = HashMap!(string, int)(256);

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
                    version(avl)
                        auto tree = this.chainsByContig.require(c.targetName, new IntervalAVLTree!ChainLink);
                    version(splay)
                        auto tree = this.chainsByContig.require(c.targetName, new IntervalSplayTree!ChainLink);
                    
                    // Insert all intervals from the chain into the tree
                    foreach(link; *c.links)
                    {
                        version(avl)    { uint cnt; (*tree).insert( new IntervalTreeNode!ChainLink(*link), cnt ); }
                        version(splay)  (*tree).insert(*link);
                    }

                    // Record the destination ("query") contig needed for VCF header
                    this.qContigSizes[c.queryName] = c.querySize;
                }
                chainStart = i;
            }
        }
        // !!! Don't forget the last chain !!!
        auto c = Chain(chainArray[chainStart .. $]);
        debug hts_log_trace(__FUNCTION__, format("Final Chain: %s", c));

        // Does this contig exist in the map?
        version(avl)
            auto tree = this.chainsByContig.require(c.targetName, new IntervalAVLTree!ChainLink);
        version(splay)
            auto tree = this.chainsByContig.require(c.targetName, new IntervalSplayTree!ChainLink);

        // Insert all intervals from the chain into the tree
        foreach(link; *c.links)
        {
            version(avl)    { uint cnt; (*tree).insert( new IntervalTreeNode!ChainLink(*link), cnt ); }
            version(splay)  (*tree).insert(*link);
        }

        // Record the destination ("query") contig needed for VCF header
        this.qContigSizes[c.queryName] = c.querySize;

        debug { // debug-only to speed startup
            foreach(contig; this.chainsByContig.byKey) {
                hts_log_trace(__FUNCTION__, format("Contig: %s", contig));
            }
        }

        // show me what's in the last accessed tree:
        debug { // debug-only to speed startup
            version(avl) {}
            else {
                while(tree.iteratorNext() !is null)
                hts_log_trace(__FUNCTION__, format("sorted tree entry: %s", *tree.cur));
            }
        }  
    }   // END ChainFile ctor

    /** Lift a single coordinate (in zero-based, half-open system),
        mutating function parameters 
        
        Returns:    number of results (0 or 1)
    */
    int liftDirectly(ref string contig, ref int coord)
    {
        auto i = BasicInterval(coord, coord + 1);
        auto o = this.chainsByContig[contig].findOverlapsWith(i);   // returns Node*

        const auto nres = o.length;
        if (!nres) return 0;
        else {
            const auto isect = o.front().interval.intersect(i);

            // interval is type ChainLink
            contig = isect.qContig;
            coord = isect.qStart;
            return 1;
        }
    }

    /** Lift coordinates from one build to another

        TODO: error handling (at cost of speed)

        Returns:    array (TODO, range) of ChainLink
                    Initially was ChainLink, then ChainInterval, then ChainLink again
                    because we need strandInvert *and* delta *and* source/dest
    */
    ChainLink[] lift(string contig, int start, int end) // can't be const method since findOverlapsWith mutates tree
    {
        auto i = BasicInterval(start, end);
        auto o = this.chainsByContig[contig].findOverlapsWith(i);   // returns Node*(s)

        // marked as debug because in hot code path
        debug foreach(x; o) hts_log_debug(__FUNCTION__, format("%s", *x));

        if (o.length == 0)  // no match
        {
            debug hts_log_debug(__FUNCTION__, "No match to interval");

            // -O3 will elide the stack alloc, thanks godbolt.org !
            ChainLink[] ret;
            return ret;
        }
        else if (o.length == 1) // one match
        {
            debug hts_log_debug(__FUNCTION__, "One match to interval");

            // intersect makes the chain link comply with bounds of interval
            const auto isect = o.front().interval.intersect(i);
            //TODO here is where we would want to report truncated interval

            /+
            debug if(isect.start + o.front().interval.delta == 248_458_169)
                hts_log_debug(__FUNCTION__, format("%s", o.front().interval));
            +/
            return [isect];
        }
        else    // TODO optimize; return Range
        {
            debug hts_log_debug(__FUNCTION__, "Multiple matches to interval");

            // [] needed for UnrolledList (opSlice to return Range over the container); was not needed when dynamic array
            auto isect = o[].map!(x => x.interval.intersect(i));

            return array(isect);
        }
    }
}

/// Return the intersection of a ChainLink with a generic interval
/// Comparison and trimming is on + strand only, but does support
/// ChainLink with query(dest) on the - strand
pragma(inline, true)
ChainLink intersect(IntervalType)(ChainLink link, const IntervalType other)
if (__traits(hasMember, IntervalType, "start") &&
    __traits(hasMember, IntervalType, "end"))
{
    // Trimming is "inwards", i.e., ===== -> -===-
    const auto trimmedStart = max(link.start, other.start);
    const auto trimmedEnd   = min(link.end,   other.end);

    const auto delta = trimmedStart - link.start;
    link.qStart += link.invert * delta;

    link.tStart = trimmedStart;
    link.tEnd   = trimmedEnd;

    return link;
}

/// return the intersection of a ChainLink with an interval
/// NOTE: this requires the first two elements of InteravalType
/// be (start, end) or have ctor(start, end)
/// TODO: error handling (at cost of speed)
pragma(inline, true)
ChainLink intersectX(IntervalType)(ChainLink int1, IntervalType int2)
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
            ChainInterval(int1.target.contig, trimmedStart, trimmedEnd, int1.target.strand),
            ChainInterval(int1.query.contig,  int1.query.start + startDiff,
                                            int1.query.end - endDiff,
                                            int1.query.strand),
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
*/

/// swap start, end if start > end
@safe
@nogc nothrow
void orderStartEnd(ref int start, ref int end)
{
    int s;
    s = start;
    start = min(start, end);
    end = max(s, end);
}