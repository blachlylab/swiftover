module liftover;

import std.algorithm : splitter;
import std.container.rbtree;
import std.container.util : make;
import std.format;

/// Represents a query or target interval
/// zero based, half-open coordinates
/// If conver to class, make implement Interface Interval (need: start, end, overlaps, opCmp)
struct ChainInterval
{
	string contig;
	uint start;
	uint end;

	/** Detect overlap between this interval and other given interval

	Only defined if intervals on same contig

	return true in any of the following four situations:
		self   =====    =======
		other =======  =======
		
		self  =======  =======
		other   ===      =======

	return false in any other scenario:
		self  =====         |        =====
		other        =====  |  =====

	*/
	bool overlaps(ref const ChainInterval other) const
	{
		// Intervals frmo different contigs are non-overlapping
		//if (this.contig != other.contig) return false;	// possibly swap for assert?
		assert(this.contig == other.contig);

		// self   =====    =======
		// other =======  =======
		if (other.start <= this.start &&  this.start <= other.end) return true;

		// self  =======  =======
		// other   ===      =======
		else if (this.start <= other.start && other.start <= this.end) return true;

		// self  =====         |        =====
		// other        =====  |  =====
		else return false;
	}

	int opCmp(ref const ChainInterval other) const
	{
		// Intervals from different contigs are incomparable
		assert(this.contig == other.contig);	// formerly return 0, but implies equality, leads to "duplicate insert" bug
		
		if (this.start < other.start) return -1;
		else if(this.start > other.start) return 1;
		else if(this.start == other.start && this.end < other.end) return -1;	// comes third as should be less common
		else if(this.start == other.start && this.end > other.end) return 1;
		else return 0;	// would be reached in case of equality (although we do not expect)
	}

	string toString() const
	{
		return format("%s:%d-%d", this.contig, this.start, this.end);
	}
}

/// Represents a mapping from ChainInterval -> ChainInterval
/// Beacuse this is what's stored in RBTree nodes,
/// should also impoement Interface Interval, but delegate to member 'query'
struct ChainBlock
{
	ChainInterval query;
	ChainInterval target;

	/// Overlap defined in terms of query interval
	bool overlaps(ref const ChainBlock other) const
	{
		return this.query.overlaps(other.query);
	}

	/// Overload <, <=, >, >=
	/// ChainBlocks should be ordered [in the interval tree] according to the query
	int opCmp(ref const ChainBlock other) const
	{
		return query.opCmp(other.query);
	}

	string toString() const
	{
		return format("%s → %s", this.query, this.target);
	}
}

/// The Chain type represents a UCSC chain object, including all
/// the fields from the header line and each block of mappings
/// for that chain.
struct Chain(R)
if(isInputRange!R)
{
	long Score;

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

	string id;

	ChainBlock[] blocks;

	this(string chainHeader, R lines)
	{
		//chain 20851231461 chr1 249250621 + 10000 249240621 chr1 248956422 + 10000 248946422 2
		// assumes no errors in chain line
		auto s = chainHeader.splitter(',');
		
	}

	string toString() const
	{
		return format("%s:%d-%d to %s:%d-%d :: %d ChainBlocks", this.queryName, this.queryStart, this.queryEnd,
			this.targetName, this.targetStart, this.targetEnd, this.blocks.length);
	}
}

struct LiftoverGroup
{
	string target;	/// destination
	string query;	/// source

	RedBlackTree!ChainInterval[string] contigs;
}


int loadChainFile(string target, string query, string fn)
{
	return 0;
}

unittest
{
	import std.stdio : writeln, writefln;

	auto ci1 = ChainInterval("chr1", 0, 10);
	auto ci2 = ChainInterval("chr1", 10, 20);

	auto ci3 = ChainInterval("chr1", 100, 200);
	auto ci4 = ChainInterval("chr1", 200, 300);

	auto cb1_2 = ChainBlock(ci1, ci2);
	auto cb3_4 = ChainBlock(ci3, ci4);

	RedBlackTree!ChainBlock rbt = make!RedBlackTree(cb1_2, cb3_4);

	writefln("Length of rbt: %d", rbt.length);

	writeln(rbt);

	foreach(e; rbt) {
		writeln(e);
	}
}
