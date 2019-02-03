module liftover;

import std.algorithm : splitter;
import std.container.rbtree;
import std.format;

/// Represents a query or target interval
/// zero based, half-open coordinates
struct ChainInterval
{
	string contig;
	uint start;
	uint end;

	int opCmp(ref const ChainInterval other) const
	{
	}

	string toString() const
	{
		return format("%s:%d-%d", this.contig, this.start, this.end);
	}
}

/// Represents a mapping from ChainInterval -> ChainInterval
struct ChainBlock
{
	ChainInterval query;
	ChainInterval target;

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

	RedBlackTree!ChainBlock[string] contigs;
}


int loadChainFile(string target, string query, string fn)
{
	return 0;
}

