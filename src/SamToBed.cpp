#include <iostream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>
#include <string.h>
#include "LiftOver.h"

using namespace std;

void PrintUsage() {
	cout << "Usage: samToBed file.sam [--reportAccuracy] [--reportDivergence] [--blocks]" << endl;
	cout << "    --reportAccuracy   Reports accuracy as #matches/(#matches+#ins-bases+#del-bases+#mismatch) " << endl
			 << "    --reportIdentity   Reports identity as #matches/(#matches+#ins-events+#del-events+#mismatch)" << endl
			 << "    --blocks           Write contiguously aligned blocks." << endl
			 << "    --minLength        Only print from sequences at least this length" << endl
			 << "    --seqLength        Print the length of the query sequence" << endl
			 << "    --ignoreN          Do not penalize \"N\"'s" << endl
			 << "    --useXS            Use XS when computing read pos." << endl
			 << "    --flag             Output the sam flag." << endl
			 << "    --useH             Incorporate hard clipping when computing the boundaries of the reads." << endl;	
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		PrintUsage();
		exit(1);
	}

	ifstream samIn(argv[1]);
	int argi = 2;
	bool reportAccuracy = false;
	bool writeBlocks = false;
	bool reportDivergence = false;
	int minLength = 0;
	bool printSeqLength = false;
	bool penalizeN = true;
	bool useXS = false;
	bool useH = false;
	bool writeFlag = false;

	while (argi < argc) {
		if (strcmp(argv[argi], "--reportAccuracy") == 0) {
			reportAccuracy = true;
		}
		else if (strcmp(argv[argi], "--reportIdentity") == 0) {
			reportDivergence = true;
		}
		else if (strcmp(argv[argi], "--blocks") == 0) {
			writeBlocks = true;
		}
		else if (strcmp(argv[argi], "--seqLength") == 0) {
			printSeqLength = true;
		}
		else if (strcmp(argv[argi], "--ignoreN") == 0) {
			penalizeN = false;
		}
		else if (strcmp(argv[argi], "--minLength") == 0) {		
			++argi;
			minLength = atoi(argv[argi]);
		}
		else if (strcmp(argv[argi], "--useXS") == 0) {
			useXS=true;
		}
		else if (strcmp(argv[argi], "--flag") == 0) {
			writeFlag=true;
		}
		else if (strcmp(argv[argi], "--useH") == 0) {
			useH=true;
		}

		else {
			PrintUsage();
			cout << "Unknown option " << argv[argi] << endl;
			exit(1);
		}
		++argi;
	}

	while (samIn) {
		string line;
		getline(samIn, line);
		
		if (line.size() > 0 and line[0] == '@') {
			continue;
		}
		if (line.size() == 0) {
			continue;
		}

		stringstream lineStrm(line);
		string contig;
		int flag;
		string chrom;
		int pos;
		int mapq;
		string cigar;
		string t1,t2;
		int alnLength;
		string seq;
		int strandFlag = 0x16;
		lineStrm >> contig >> flag >> chrom >> pos >> mapq >> cigar >> t1 >> t2 >> alnLength >> seq;
		int strand = 0;
		int nStrand = 1;
		if (flag & strandFlag) {
			strand = 1;
			nStrand=-1;
		}
		if (chrom == "*") {
			continue;
		}
		if (seq.size() < minLength) {
			continue;
		}
		string cigarPrefix = cigar.substr(0,100);
		string cigarSuffix = cigar.substr(max(0, (int)cigar.size()-100), 100);
		int s1 = 0;
		int nMatch = 0, nIns = 0, nDel = 0, nMisMatch = 0;
		int nInsEvent = 0, nDelEvent = 0;
		int frontClip = 0, endClip = 0;
		int i = 0;
		bool alnStarted = false;
		int refPos = pos;
		int qPos = 0;
		int xs=0,xe=0,xl=0;
		if (useXS) {
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
		if (useXS) {
			if (nStrand == 1) {
				frontClip = xs-1;
			}
			else {
				frontClip = xl - xe;
			}
		}
		
		while (i < cigar.size()) {
			int len = atoi(&cigar[i]);
			while (i < cigar.size() and cigar[i] >= '0' and cigar[i] <= '9') { i++;}
			char op = cigar[i];
			 if (op == 'S' or op == 'H') {
				 if (op == 'H' and useH and alnStarted == false) {
					 frontClip += len;
				 }
				 if (op == 'S' and alnStarted == false and useXS == false) {					 
					frontClip += len;					
					qPos += len;					
				}
				else if ((op == 'H' or op == 'S') and alnStarted == true) {
					endClip = len;
				}
			}
			else {
				alnStarted = true;
				if ((op == 'M' or op == 'X' or op == '=') and writeBlocks)  {
					cout << chrom << "\t" << refPos << "\t" << refPos + len << endl;
				}
					
				if (op == 'M') {
					nMatch += len;
					refPos += len;
				}
				if (op == '=') {
					nMatch +=len;
					refPos +=len;
				}
				else if (op == 'X') {
					
					nMisMatch+=len;
					refPos+=len;
					int p = qPos;
					int numN = 0;
					for (p = qPos; p < qPos + len; p++) {
						if (seq[p] == 'N') numN ++;
					}
					if (penalizeN == false) {
						nMisMatch -= numN;
					}
				}
				else if (op == 'I') {
					nIns += len;
					nInsEvent++;
				}
				else if (op == 'D') {
					nDel += len;
					refPos +=len;
					nDelEvent++;
				}
			}
			i+=1;
		}
		int denom = nMatch + nMisMatch + nIns + nDel;
		float identity = 0;
		if (denom > 0) {
			if (reportDivergence == false) {
				identity = ((float) nMatch) / denom;
			}
			else {
				identity = ((float)nMatch)/(nMatch + nInsEvent + nDelEvent + nMisMatch);			
			}
		}
		if (writeBlocks == false) {
			cout << chrom << "\t" << pos - 1 << "\t" << pos -1 + nMatch + nMisMatch + nDel  << "\t" 
					 << contig << "\t" << strand << "\t" << frontClip << "\t" << frontClip + nMatch + nMisMatch + nIns << "\t" 
					 << mapq << "\t";
			if (reportAccuracy or reportDivergence) {

				cout << identity << "\t"  
						 << nMatch << "\t" << nMisMatch << "\t"
						 << nIns << "\t" << nDel;
			}
			if (printSeqLength) {
				cout << "\t" << seq.size();
			}
			if (writeFlag) {
				cout << "\t" << flag;
			}
			cout << endl;
		}
	}
}
	
