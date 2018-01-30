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

	string s;
	vector<pair<Hit, string>> hits;
	while (getline(fin, s)) {
		string cigar;
		Hit h = Hit::from_bed(s, &cigar);
		assert(h.ref != nullptr);
		assert(h.query != nullptr);
		if (tie(h.query->name, h.query_start, h.query_end) > tie(h.ref->name, h.ref_start, h.ref_end)) {
			swap(h.query->name, h.ref->name);
			swap(h.query_start, h.ref_start);
			swap(h.query_end, h.ref_end);
		}
		hits.push_back({h, cigar});
	}
	sort(hits.begin(), hits.end(), [](const pair<Hit, string> &a, const pair<Hit, string> &b) {
		return 
			tie(a.first.ref->is_rc, a.first.query->name, a.first.ref->name, a.first.query_start, a.first.ref_start) <
			tie(b.first.ref->is_rc, b.first.query->name, b.first.ref->name, b.first.query_start, b.first.ref_start);
	});

	int hit_count = 0, out_count = 0;
	string prev;
	for (auto &hs: hits) {
		auto &h = hs.first;
		auto &cigar = hs.second;
		string fa = fr.get_sequence(h.query->name, h.query_start, &h.query_end);
		string fb = fr.get_sequence(h.ref->name, h.ref_start, &h.ref_end);
		assert(!h.query->is_rc); 
		if (h.query->is_rc) {
			fa = rc(fa);
		}
		if (h.ref->is_rc) {
			fb = rc(fb);
		}
		assert(cigar.size());
		h.aln = Alignment(fa, fb, cigar);
		// h.aln.trim();
		// hits.push_back(h);
	// }
	// eprn("Loaded {} hits", hits.size());
	// // hits = merge(hits, 250);
	// // eprn("After merging remaining {} hits", hits.size());
	// // for (auto &h: hits) {
	// // 	prn("{}", h.to_bed(false));
	// // }
	// // exit(0);
	// for (auto &h: hits) {
		int align_length = h.aln.span();
		int indel_a = 0;
		int indel_b = 0;
		int alignB = 0;
		int matchB = 0;
		int mismatchB = 0;
		int transitionsB = 0;
		int transversionsB = 0;

		int uppercaseA = 0;
		int uppercaseB = 0;

		for (int i = 0; i < align_length; i++) {
			char a = toupper(h.aln.align_a[i]);
			char b = toupper(h.aln.align_b[i]);
			indel_a += a == '-';
			indel_b += b == '-';
			matchB += a != '-' && a == b;
			uppercaseA += (h.aln.align_a[i] != '-' && h.aln.align_a[i] != 'N' && isupper(h.aln.align_a[i]));
			uppercaseB += (h.aln.align_b[i] != '-' && h.aln.align_b[i] != 'N' && isupper(h.aln.align_b[i]));
			if (a != '-' && b != '-') {
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
		
		// TODO handle this smarter
		bool same_chr = h.query->name == h.ref->name && h.query->is_rc == h.ref->is_rc;
		int overlap = !same_chr ? 0 : max(0, 
			min(h.query_end, h.ref_end) - max(h.query_start, h.ref_start));
		bool too_big_overlap = 
			(h.query_end - h.query_start - overlap) < 100 ||
			(h.ref_end - h.ref_start - overlap) < 100;
		too_big_overlap &= same_chr;
		// too_big_overlap = 0;

		double errorScaled = (h.aln.gaps() + h.aln.mismatches()) / 
			double(h.aln.gaps() + h.aln.mismatches() + h.aln.matches());

		if (uppercaseA >= 100 && uppercaseB >= 100 && !too_big_overlap) {
			string l = h.to_bed(false);
			
			if (l != prev) {
				h.name = fmt::format("S{:05}", ++out_count);
				prn("{}\t"
					"{}\t{}\t{}\t{}\t{}\t{}\t"
					"{}\t{}\t{}\t{}\t{}\t{}\t"
					"{}\t{}\t{}\t"
					"{}\t{}\t{}\t{}", 
					h.to_bed(false), 
					align_length, indel_a, indel_b, alignB, matchB, mismatchB,
					transitionsB, transversionsB, fracMatch, fracMatchIndel, jcK, k2K,
					h.aln.gaps(),
					uppercaseA, uppercaseB,
					h.aln.matches(), h.aln.mismatches(), h.aln.gaps(), h.aln.gap_bases()
				);
				if (out_count > 10) exit(0);
			}
			prev = l;
		}
		eprnn("\rProcessed hit {:n}", ++hit_count);
		// exit(0);
	}
	eprn("\nRead {:n} hits, wrote {:n} SDs", hit_count, out_count);
	eprn("\nDone!");
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
