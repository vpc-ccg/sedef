/// 786 

/******************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>

#include "align.h"
#include "common.h"
#include "fasta.h"
#include "hit.h"
#include "merge.h"

using namespace std;

/******************************************************************************/

void stats(const string &ref_path, const string &bed_path) 
{
	FastaReference fr(ref_path);
	ifstream fin(bed_path.c_str());
	if (!fin.is_open()) {
		throw fmt::format("BED file {} does not exist", bed_path);
	}

	string line;
	vector<Hit> hits;
	while (getline(fin, line)) {
		string cigar;
		auto h = Hit::from_bed(line, &cigar);
		// string fa = fr.get_sequence(h.query->name, h.query_start, h.query_end);
		// string fb = fr.get_sequence(h.ref->name, h.ref_start, h.ref_end);
		// assert(!h.query->is_rc);
		// if (h.query->is_rc) {
		// 	fa = rc(fa);
		// }
		// if (h.ref->is_rc) {
		// 	fb = rc(fb);
		// }
		// assert(cigar.size());
		// h.aln = align(fa, fb, cigar);
		hits.push_back(h);
	}
	eprn("Loaded {} hits", hits.size());
	hits = merge(hits, 250);
	eprn("After merging remaining {} hits", hits.size());
	for (auto &h: hits) {
		prn("{}", h.to_bed(false));
	}
	exit(0);
	for (auto &h: hits) {
		int align_length = h.aln.span();
		int indel_a = 0;
		int indel_b = 0;
		int alignB = 0;
		int matchB = 0;
		int mismatchB = 0;
		int transitionsB = 0;
		int transversionsB = 0;

		int lowercaseA = 0;
		int lowercaseB = 0;

		for (int i = 0; i < align_length; i++) {
			char a = toupper(h.aln.align_a[i]);
			char b = toupper(h.aln.align_b[i]);
			indel_a += a == '-';
			indel_b += b == '-';
			matchB += a != '-' && a == b;
			if (a != '-' && b != '-') {
				lowercaseA += islower(a);
				lowercaseB += islower(b);
				alignB += 1;
				if (a != b) {
					mismatchB += 1;
					if (a == 'A' || a == 'G' ) {
						transitionsB += b == 'A' || b == 'G';
						transversionsB += !(b == 'A' || b == 'G');
					} else {
						transitionsB += b == 'C' || b == 'T';
						transversionsB += !(b == 'C' || b == 'T');
					}
				}
			}
		}

		double fracMatch = double(matchB) / (alignB);
		double fracMatchIndel = double(matchB) / (align_length);
	
		double jcp = double(mismatchB) / (alignB);
		double jcK = -0.75 * log(1.0 - 4.0 / 3 * jcp);
		
		double p = double(transitionsB) / (alignB);
		double q = double(transversionsB) / (alignB);
		double w1 = 1.0 / (1 - 2.0 * p - q);
		double w2 = 1.0 / (1 - 2.0 * q);
		double k2K = 0.5 * log(w1) + 0.25 * log(w2);
		
		prn("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
			h.to_bed(false), 
			align_length, indel_a, indel_b, alignB, matchB, mismatchB,
			transitionsB, transversionsB, fracMatch, fracMatchIndel, jcK, k2K,
			lowercaseA, lowercaseB);
	}
}

/******************************************************************************/

void stats_main(int argc, char **argv) 
{
	if (argc < 2) {
		eprn("invalid usage");
		return;
	}
	string ref_path = argv[0];
	string bed_path = argv[1];
	stats(ref_path, bed_path);
}
