/// 786

/******************************************************************************/

#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>

#include "align.h"
#include "common.h"
#include "fasta.h"
#include "extern/ksw2.h"

using namespace std;

/******************************************************************************/

alignment_t align_helper(const string &tseq, const string &qseq, int sc_mch, int sc_mis, int gapo, int gape)
{
	const int STEP = 50 * 1000; // Max. alignment size (if larger, split into pieces)

	int8_t a = (int8_t)sc_mch, b = sc_mis < 0 ? (int8_t)sc_mis : (int8_t)(-sc_mis); // a>0 and b<0
	int8_t mat[25] = { 
		a, b, b, b, 0, 
		b, a, b, b, 0, 
		b, b, a, b, 0, 
		b, b, b, a, 0, 
		0, 0, 0, 0, 0 
	};
	deque<pair<char, int>> cigar;
	for (int SP = 0; SP < min(tseq.size(), qseq.size()); SP += STEP) {
		ksw_extz_t ez;
		ksw_extz2_sse(0, 
			min(STEP, (int)(qseq.size() - SP)), (const uint8_t*)(qseq.c_str() + SP), 
			min(STEP, (int)(tseq.size() - SP)), (const uint8_t*)(tseq.c_str() + SP), 
			5, mat, // M; MxM matrix
			gapo, gape, 
			-1, -1, // band width; off-diagonal drop-off to stop extension (-1 to disable)
			0, &ez);
		for (int i = 0; i < ez.n_cigar; ++i) 
			cigar.push_back(make_pair("MID"[ez.cigar[i]&0xf], ez.cigar[i]>>4));
		free(ez.cigar);
	}
	return alignment_t{ 0, 0, qseq.size(), 0, 0, tseq.size(), qseq, tseq, "", "", "", cigar, 0 };
}

alignment_t align(string fa, string fb, int match, int mismatch, int gap_open, int gap_extend) 
{
	string xa = fa, xb = fb;
	for (auto &c: fa) c = align_dna(c);
	for (auto &c: fb) c = align_dna(c);
	auto a = align_helper(fa, fb, match, mismatch, gap_open, gap_extend);
	a.a = xa, a.b = xb;
	return a;
}

/******************************************************************************/

void align_ref(FastaReference &fr, int line, const string &s, chrono::time_point<chrono::high_resolution_clock> &t_start)
{
	vector<string> ss = split(s, '\t');

	char ca[500]; int sa, ea;
	char cb[500]; int sb, eb;
	char sta[10], stb[10];
	sscanf(s.c_str(), "%s %d %d %s %d %d %*s %s %s", ca, &sa, &ea, cb, &sb, &eb, sta, stb);

	// sta[0] = '+';
	// sscanf(s.c_str(), "%s %d %d %*s %*s %s %s %d %d", ca, &sa, &ea, stb,  cb, &sb, &eb);
	// stb[0] = stb[0] == '_' ? '-' : '+';
 
	assert(sta[0] == '+');
	bool rc = (stb[0] == '-');

	string fa = fr.getSubSequence(ca, sa, ea - sa);
	string fb = fr.getSubSequence(cb, sb, eb - sb);
	if (rc) {
		reverse(fb.begin(), fb.end());
		transform(fb.begin(), fb.end(), fb.begin(), rev_dna);
	}

	auto aln = align(fa, fb);

	eprn("\t{}s\t{}", 
		chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t_start).count() / 1000.00, 
		line);
	cout << s << "\t" << aln.cigar_string() << "\t" << fa << "\t" << fb << endl;
}

void align_main(string ref_path, string bed_path, int resume_after) 
{
	auto start = chrono::high_resolution_clock::now();
	
	FastaReference fr(ref_path);
	ifstream fin(bed_path.c_str());

	int nl = 0;
	string s;
	while (getline(fin, s)) {
		if (nl > resume_after) align_ref(fr, nl, s, start);
		nl++;
	}

	eprn("Finished bucket {}: read {} lines", string(bed_path), nl);
	eprnn("Done!\n");
}

/******************************************************************************/

string alignment_t::cigar_string()
{
	string res;
	for (auto &p: cigar) res += fmt::format("{}{}", p.second, p.first);
	return res;
}

alignment_t alignment_t::trim() 
{
	auto result = *this;
	while (result.cigar.size()) {
		if (result.cigar[0].first == 'D') {
			result.a = result.a.substr(result.cigar[0].second);
			result.start_a += result.cigar[0].second;
			result.cigar.pop_front();
			continue;
		}
		if (result.cigar[0].first == 'I') {
			result.b = result.b.substr(result.cigar[0].second);
			result.start_b += result.cigar[0].second;
			result.cigar.pop_front();
			continue;
		}
		if (result.cigar.back().first == 'D') {
			result.end_a -= result.cigar.back().second;
			result.a = result.a.substr(0, result.a.size() - result.cigar.back().second);
			result.cigar.pop_back();
			continue;
		}
		if (result.cigar.back().first == 'I') {
			result.end_b -= result.cigar.back().second;
			result.b = result.b.substr(0, result.b.size() - result.cigar.back().second);
			result.cigar.pop_back();
			continue;
		}
		break;
	}

	result.populate_nice_alignment();
	return result;
}

void alignment_t::populate_nice_alignment()
{
	align_a = "";
	align_b = "";
	alignment = "";
	int ia = 0, ib = 0;
	for (auto &c: cigar) {
		for (int i = 0; i < c.second; i++) {
			assert(c.first != 'M' || ia < a.size());
			assert(c.first != 'M' || ib < b.size());
			if (c.first == 'M' && toupper(a[ia]) == toupper(b[ib])) {
				alignment += "|";
			} else {
				alignment += "*";
			}
			if (c.first != 'D') align_b += b[ib++];
			else                align_b += "-";
			if (c.first != 'I') align_a += a[ia++];
			else                align_a += "-";
		}
	}
}

string alignment_t::print(int width)
{
	if (!alignment.size()) populate_nice_alignment();

	string res;
	int qa = start_a;
	int qb = start_b;


	int gaps = 0, gap_bases = 0, mismatches = 0, matches = 0;
	for (auto &c: cigar) {
		if (c.first != 'M') gaps++, gap_bases += c.second;
	}
	for (int i = 0; i < alignment.size(); i += width) {
		if (align_a[i] != '-' && align_b[i] != '-') {
			if (toupper(align_a[i]) == toupper(align_b[i])) {
				matches++;
			} else {
				mismatches++;
			}
		}
	}

	res += fmt::format(
		"       A: {:>9}..{:<9} (len {:7})    Gaps:       {:5} = {:.0f}% ({})\n"
		"       B: {:>9}..{:<9} (len {:7})    Mismatches: {:5} = {:.0f}%\n"
		"   CIGAR: {}\n",
		start_a, end_a, end_a - start_a, 
		gap_bases, 100.0 * gap_bases / double(alignment.size()), gaps,
		start_b, end_b, end_b - start_b,
		mismatches, 100.0 * mismatches / double(alignment.size()),
		cigar_string()
	);
	for (int i = 0; i < alignment.size(); i += width) {
		res += fmt::format(
			"   {:10}: {} {}\n   {:10}  {}\n   {:10}: {}\n", 
			qa, align_a.substr(i, width), i+align_a.substr(i, width).size(),
			"", alignment.substr(i, width), 
			qb, align_b.substr(i, width));
		for (auto c: align_a.substr(i, width)) if (c != '-') qa++;
		for (auto c: align_b.substr(i, width)) if (c != '-') qb++;
		res += "\n";
	}
	return res;
}

void alignment_t::cigar_from_alignment() 
{
	cigar.clear();
	int sz = 0;
	char op = 0, top;
	for (int i = 0; i < alignment.size(); i++) {
		if (align_a[i] == '-') {
			top = 'I';
		} else if (align_b[i] == '-') {
			top = 'D';
		} else {
			top = 'M';
		}

		if (op != top) {
			if (op) cigar.push_back(make_pair(op, sz));
			op = top, sz = 0;
		}
		sz++;
	}
	cigar.push_back(make_pair(op, sz));
}

vector<alignment_t> alignment_t::max_sum(int min_span)
{
	assert(alignment.size());

	// start/end in the alignment
	int best_score = -1, best_start = 0, best_end = 0;
	int score_so_far = 0, start_so_far = 0, end_so_far = 0;

	vector<int> lookup_a(alignment.size() + 1, 0); int ia = 0;
	vector<int> lookup_b(alignment.size() + 1, 0); int ib = 0;

	vector<pair<int, int>> hits; // start, end
	for (int i = 0; i < alignment.size(); i++) {
		int score = (alignment[i] == '|' ? 1 : -1);

		if (score > score_so_far) {
			if (best_end - best_start >= min_span 
				&& (!hits.size() || hits.back() != make_pair(best_start, best_end))) 
			{
				hits.push_back({best_start, best_end});
			}
			score_so_far = score;
			start_so_far = i;
			end_so_far = i + 1;
			best_score = score_so_far;
			best_start = start_so_far, best_end = end_so_far;
		} else {
			score_so_far += score;
			end_so_far++;
		}

		if (score_so_far > best_score) {
			best_score = score_so_far;
			best_start = start_so_far, best_end = end_so_far; 
		}

		lookup_a[i] = ia;
		lookup_b[i] = ib;
		if (align_a[i] != '-') ia++; 
		if (align_b[i] != '-') ib++;
	}
	lookup_a[alignment.size()] = ia;
	lookup_b[alignment.size()] = ib;

	if (best_end - best_start >= min_span 
		&& (!hits.size() || hits.back() != make_pair(best_start, best_end))) 
	{
		hits.push_back({best_start, best_end});
	}

	vector<alignment_t> results;
	for (auto &h: hits) {
		// translate: a --> b
		results.push_back({
			chr_a, lookup_a[h.first], lookup_a[h.second],
			chr_b, lookup_b[h.first], lookup_b[h.second],
			a.substr(lookup_a[h.first], lookup_a[h.second] - lookup_a[h.first]),
			b.substr(lookup_b[h.first], lookup_b[h.second] - lookup_b[h.first]),
			align_a.substr(h.first, h.second - h.first),
			align_b.substr(h.first, h.second - h.first),
			alignment.substr(h.first, h.second - h.first),
			{}, 0
		});
		results.back().cigar_from_alignment();
	}

	return results;
}
