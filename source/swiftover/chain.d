module swiftover.chain;

import std.array : appender;
import std.algorithm : splitter;
import std.ascii : isWhite;
import std.format;

/// Represents a mapping from ChainInterval -> ChainInterval
/// Previously, this contained two "ChainInterval" structs (one target, one query)
/// But they are merged for efficiency
/// Beacuse this is what's stored in RBTree nodes,
/// should also implement Interface Interval, but delegate to member 'query'
/// (tree sorted by query, not target)
struct ChainLink
{
    int qcid;   /// Query contig id
    int qstart; /// Query start
    int qend;   /// Query end

    int tcid;   /// Target contig id
    int tstart; /// Target start
    int tend;   /// Target end

    // To work with the interval tree, and overlap functions
    alias start = qstart;
    alias end = qend;

    /// Overload <, <=, >, >= for ChainLink/ChainLin; compare query
	int opCmp(ref const ChainLink other) const
	{
		// Intervals from different contigs are incomparable
		assert(this.qcid == other.qcid);	// formerly return 0, but implies equality, leads to "duplicate insert" bug when using std.container.rbtree
		
		if (this.start < other.start) return -1;
		else if(this.start > other.start) return 1;
		else if(this.start == other.start && this.end < other.end) return -1;	// comes third as should be less common
		else if(this.start == other.start && this.end > other.end) return 1;
		else return 0;	// would be reached in case of equality (although we do not expect)
	}
	/// Overload <, <=, >, >= for ChainLink/int 
	int opCmp(const int x) const
	{
		if (this.start < x) return -1;
		else if (this.start > x) return 1;
		else return 0;
	}

    // Todo, figure out if we can look up the contig id -> string mapping (static member map?)
	string toString() const
	{
		return format("(contig#%d):%d-%d → (contig#%d):%d-%d",
                    this.qcid, this.qstart, this.qend,
                    this.tcid, this.tstart, this.tend);
	}
}


/// The Chain type represents a UCSC chain object, including all
/// the fields from the header line and each block of mappings
/// for that chain.
/// 
/// Specifications: https://genome.ucsc.edu/goldenpath/help/chain.html
struct Chain(R)
if(isInputRange!R)
{
	long score;         /// Chain score

	string targetName;
	long targetSize;
	char targetStrand;	/// +,-
	long targetStart;
	long targetEnd;
	
	string queryName;
	long querySize;
	char queryStrand;
	long queryStart;
	long queryEnd;

	int id;             /// chain id

	ChainLink*[] links;  /// query and target intervals in 1:1 bijective relationship

    static auto hfields = appender!(char[][]);  /// header fields; statically allocated to save GC allocations
    static auto dfields = appender!(int);       /// alignment data fields; statically allocated to save GC allocations

    /// Construct Chain object from a header and range of lines comprising the links or blocks
	this(string chainHeader, R lines)
	{
        // Example chain header line: 
		// chain 20851231461 chr1 249250621 + 10000 249240621 chr1 248956422 + 10000 248946422 2
		// assumes no errors in chain line
        this.fields.clear;
        this.fields.put(chainHeader.splitter(" "));
        assert(this.fields[0] == "chain");
        this.score = this.fields[1].to!long;
        
        this.targetName = this.fields[2];
        this.targetSize = this.fields[3].to!long;
        this.targetStrand = this.fields[4][0];
        this.targetStart = this.fields[5].to!long;
        this.targetEnd = this.fields[6].to!long;

        this.queryName = this.fields[7];
        this.querySize = this.fields[8].to!long;
        this.queryStrand = this.fields[9][0];
        this.queryStart = this.fields[10].to!long;
        this.queryEnd = this.fields[11].to!long;

        this.id = this.fields[12].to!int;


        /*  Alignment data lines:
            size dt dq
            size -- the size of the ungapped alignment
            dt -- the difference between the end of this block and the beginning of the next block (reference sequence)
            dq -- the difference between the end of this block and the beginning of the next block (query sequence)
            NOTE: The last line of the alignment section contains only one number: the ungapped alignment size of the last block.
        */
        int tFrom = this.targetStart;   // accum current coordinate, target
        int qFrom = this.queryStart;    // accum current coordinate, query
        bool done;
        foreach(line; lines)
        {
            this.dfields.clear();
            this.dfields.put(line.splitter!isWhite().map(x => x.to!int));   // TODO: benchmark splitter() 

            // set up ChainLink from alignement data line
            ChainLink* link = new ChainLink;

            link.qcid = 99; // TODO, get int id for string this.queryName
            link.qstart = qFrom;
            link.qend = qFrom + size;

            link.tid = 101; // TODO, get int id for string this.targetName
            link.tstart = tFrom;
            link.tend = tFrom + size;

            if(this.dfields.length == 1)    // last block in chain
                done = true;
            else if(this.dfields.length == 3)
            {
                assert(done is false, "Malformed alignment data blocks");
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
		return format("%s:%d-%d → %s:%d-%d :: %d links", this.queryName, this.queryStart, this.queryEnd,
			this.targetName, this.targetStart, this.targetEnd, this.links.length);
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

/// Representation of UCSC-format liftover chain file, which contains multiple alignment/liftover chains
struct ChainFile
{
    string sourceBuild; /// original assembly, e.g. hg19 in an hg19->GRCh38 liftover
    alias queryBuild = sourceBuild;
    string targetBuild; /// target assembly, e.g. GRCh38 in an hg19->GRCh38 liftover

    /// Parse UCSC-format chain file into liftover trees (one tree per source/query contig)
    this(string fn)
    {
        assert(0, "WIP resume here");
    }
}