#ifndef LIFT_OVER_H_
#define LIFT_OVER_H_
#define Q  0
#define T 1
#include <vector>
#include <map>
#include <assert.h>
#include <algorithm>
using namespace std;
class Block {
public:
	int t,q,l;
	static int dir;
	Block(int &tp, int &qp, int &lp) {
		t = tp; q = qp; l=lp;
	}
	int Pos() const {
		if (dir == 0) {
			return q;
		}
		else {
			return t;
		}
	}
	int Map() {
		if (dir == 0) {
			return t;
		}
		else {
			return q;
		}
	}
	int operator<(const int p) const {
		return Pos() < p;
	}
};


typedef vector<Block> Blocks;
typedef vector<Blocks> MultiBlocks;
typedef vector<int> Strands;
typedef map<string, MultiBlocks> PosMap;
typedef map<string, Strands> StrandMap;
typedef map<string, vector<int> > LengthMap;
typedef map<string, vector<string> > ChromMap;
typedef map<string, vector<int  > > ClipMap;

bool SetPos(const int nStrand, Blocks &blocks, int index, int pos, int &mapPos) {
	assert(blocks[index].Pos() <= pos);


	if (index == blocks.size() - 1 and pos > blocks[index].Pos() + blocks[index].l) {
		return false;
	}

	if (index < blocks.size() - 1) {
       assert(blocks[index+1].Pos() > pos);
	}

	if (blocks[index].Pos() + blocks[index].l > pos) {
       if ( nStrand ) {

          mapPos = blocks[index].Map() + blocks[index].Pos() - pos;

       }
       else {

          mapPos = blocks[index].Map() + pos - blocks[index].Pos()  ;
       }
	}
	else {
		// 
		// This position is in a gap.
       if ( nStrand ) {
          mapPos = blocks[index].Map() - blocks[index].l;
       }
       else {
          mapPos = blocks[index].Map() + blocks[index].l;
       }
	}
	return true;
}

bool SearchBlocks(const int nStrand, Blocks &blocks, int pos, int &mapPos) {
	if (blocks.size() == 0) {
		return false;
	}
	Blocks::iterator lb;
	lb = lower_bound(blocks.begin(), blocks.end(), pos);
	int i = lb - blocks.begin();
	//
    // Find the index that contains the block, or none if no such block exists.
	//

    if ( i >= blocks.size() ) { 

       // case in which lower_bound failed (all blocks are less than
       // pos).  However, the last block might still be the right now.

       i = blocks.size() - 1;

       if ( (blocks[i].Pos() + blocks[i].l) <= pos ) {
          // unfortunately, pos is greater than the last block
          return false;
       }
    }

    if ( pos < blocks[i].Pos()  ) {
       if ( i == 0 ) {
          return false;
       }
       else {
          --i;
       }
    }

    // normal case:
    if ( not ( blocks[i].Pos() <= pos and (pos < (blocks[i].Pos() + blocks[i].l) ) ) ) {
       // case in which pos is between 2 blocks due to there being
       // a deletion in this domain

       assert( i < ( blocks.size() - 1 ) );
       assert( (blocks[i].Pos() + blocks[i].l ) <= pos );
       assert( pos < blocks[i+1].Pos() );
       mapPos = blocks[i+1].Map();
       return true;
    }

		bool res = SetPos(nStrand, blocks, i, pos, mapPos);
		return res;
}

bool SearchContig(PosMap &posMap, ChromMap &chromMap, StrandMap &strandMap, LengthMap &lengthMap,
									string contig, int pos,  string &mapChrom, int &mapPos, int &mapStrand, string  &mapContig, int &mapContigIndex,
									int startBlockIndex=0) {
	
	if (posMap.find(contig) == posMap.end()) {
		return false;
	}
	int i;
	for (i = startBlockIndex; i < posMap[contig].size(); i++) {
		int searchPos = pos;

		if (SearchBlocks(strandMap[contig][i], posMap[contig][i], searchPos, mapPos) == true) {
			mapChrom = chromMap[contig][i];
			mapStrand = strandMap[contig][i];
			mapContig = contig;
			mapContigIndex = i;
			return true;
		}
	}
	return false;
}

class MapDBOptions {
 public:
	bool keepForward;
	bool useXS;
	MapDBOptions() {
		keepForward=false;
		useXS = false;
	}
};

bool GetKeyValue(string &key, string &kvp, int &value) {
	if (key.size() <= kvp.size() and kvp.substr(0,key.size()) == key) {
		string valueStr = kvp.substr(key.size()+3);
		value = atoi(valueStr.c_str());
		return true;
	}
	else {
		return false;
	}
}

bool SearchKeyValue(string key, vector<string> &kvps, int &value) {
	int i;
	for (i = 0; i< kvps.size(); i++) {
		if (GetKeyValue(key, kvps[i], value)) {
			return true;
		}
	}
	return false;
}
			
int BuildMapDB(ifstream &samIn, int dir,
							 map<string, MultiBlocks > &posMap,
							 map<string, Strands> &strands,
							 LengthMap &lengths,
							 map<string, vector<string> > &chromMap,
							 map<string, vector<string> > &seqMap,
							 ClipMap &clipMap,
							 MapDBOptions &opts) {
	int contigIndex = 0;
	string previousContig = "";
	while (samIn) {
		string line;
		getline(samIn, line);
		int tPos = 0;
		int qPos = 1; // sam files are 1-based

		if (line.size() > 0 and line[0] == '@') {
			continue;
		}
		if (line.size() == 0) {
			continue;
		}
		stringstream lineStrm(line);
		string contig;
		string flagStr;
		string chrom;
		int pos;
		int mapq;
		string cigar;
		string t1,t2;
		int alnLength;
		string seq;
		int flag;		
		lineStrm >> contig >> flag >> chrom >> pos >> mapq >> cigar >> t1 >> t2 >> alnLength >> seq;
		int xs = 0;
		int xe = 0;
		int xl = 0;
		if (opts.useXS) {
			string qual;
			lineStrm >> qual;
			vector<string> kvps;
			while(lineStrm) {
				string kvp;
				lineStrm >> kvp;
				if (kvp != "") {
					kvps.push_back(kvp);
				}
			}
			SearchKeyValue("XS", kvps, xs);
			SearchKeyValue("XE", kvps, xe);
			SearchKeyValue("XL", kvps, xl);
		}
		
		if (chrom == "*") {
			continue;
		}
		int s1 = 0;
		int nMatch = 0, nIns = 0, nDel = 0, nMisMatch = 0;
		int frontClip = 0, endClip = 0;
		int i = 0;
		bool alnStarted = false;

		string fromChrom = contig;
		string toChrom   = chrom;
		if (dir == T) {
			fromChrom = chrom;
			toChrom   = contig;
		}
		tPos = pos-1;
		int slash = fromChrom.find('/');
		if (slash != fromChrom.npos) {
		  fromChrom = fromChrom.substr(0,slash);
		}

		if (posMap.find(fromChrom) == posMap.end()) {

			posMap[fromChrom] = MultiBlocks();
			lengths[fromChrom] = vector<int>();
			chromMap[fromChrom] = vector<string>();
		}
		lengths[fromChrom].push_back(seq.size());
		posMap[fromChrom].push_back(Blocks());

		strands[fromChrom].push_back(flag & 0x10);
		
		chromMap[fromChrom].push_back(toChrom);
		seqMap[fromChrom].push_back(seq);
		Blocks *curBlocks = &(posMap[fromChrom][posMap[fromChrom].size()-1]);

		// in the case of bottom strand alignments, calculate the length
		// of the query sequence.  This requires soft clipping so the 
		// cigar will reflect the entire sequence--not just the aligned
		// portion
        
		int nStrand = 1;
		int strand = flag & 0x10;
		if ( opts.keepForward == false and flag & 0x10 ) {
			// bottom strand alignment
				nStrand = -1;
			if (opts.useXS == true) {
				qPos = xl;
			}
			else {
				int nQueryLength = 0;
				i = 0;
				while( i < cigar.size() ) {
					int len = atoi(&cigar[i]);

					// skip over the leading number
					while (i < cigar.size() and cigar[i] >= '0' and cigar[i] <= '9') { i++;}

					char op = cigar[i];


					if ( (char*) strchr( "SHMI=X", op ) != 0 ) {
						nQueryLength += len;

					}
					++i;
				}
				qPos = nQueryLength;  // start at last base position (1-based)
			}
		}
		else {
			qPos = 1; // start at 1st base position (1-based)
		}


		// end calculating length

		i = 0;  // go through cigar again

		int frontSoftClip;
		int frontHardClip = 0;
		while (i < cigar.size()) {
			int len = atoi(&cigar[i]);
			if (i == 0 and opts.useXS) {
				if (nStrand == 1) {
					frontClip = xs-1;
				}
				else {
					frontClip = xl - xe;
				}
				qPos += nStrand*frontClip;
			}
			
			while (i < cigar.size() and cigar[i] >= '0' and cigar[i] <= '9') { i++;}
			char op = cigar[i];
			
			if ((op == 'S' or op == 'H') and opts.useXS == false ) {
				if ((op == 'H' or op == 'S') and alnStarted == false) {
					if (op == 'H') {
						frontHardClip = len;
					}
					frontClip = len;
					qPos += nStrand*frontClip;
				}
					else if ((op == 'H' or op == 'S') and alnStarted == true) {
						endClip = len;
						qPos += nStrand*len;
					}
			}
			else {
				alnStarted = true;
				if (op == 'M' or op == '=' or op == 'X') {
					
					// Block.Pos() should be the *least* value of the
					// alignment in the Block. This occurs naturally
					// when the source is the reference (whether the
					// alignment is top or bottom strand), but if the
					// source is query, the the least value occurs at
					// the end of the = alignment. DG. Jan 5, 2017

					int nSetQPos = qPos;
					int nSetTPos = tPos;
					if ( ( dir == Q ) && ( flag & 0x10 ) ) {
						nSetQPos = qPos - len + 1;
						nSetTPos = tPos + len - 1;
					}

					curBlocks->push_back(Block(nSetTPos, nSetQPos, len));
					nMatch += len;
					tPos += len;
					qPos += nStrand*len;
				}
				else if (op == 'I') {
					nIns += len;
					qPos += nStrand*len;
				}
				else if (op == 'D') {
					nDel += len;
					tPos += len;
				}
			}
			i+=1;
		}

		// the Blocks must be in order according to the position they
		// will be searched by.  If they are searched by target (dir =
		// 1), they will always be in order.  If they are searched by
		// query (dir = 0), they will be in order for top strand
		// alignments, but in reverse order for bottom strand
		// alignments.  So reverse them. DG, 1/7/2017

		if ( ( dir == Q ) && ( flag & 0x10 ) ) {
			std::reverse( curBlocks->begin(), curBlocks->end() );
		}
		clipMap[fromChrom].push_back(frontHardClip);

		previousContig = contig;
		contigIndex ++;
	}
	return contigIndex;
}
int Block::dir = 0;
#endif
