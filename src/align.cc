/// 786

/******************************************************************************/

#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>
#include <glob.h>

#include "align.h"
#include "common.h"
#include "fasta.h"
#include "hit.h"
#include "extern/ksw2.h"

using namespace std;

/******************************************************************************/

Alignment Alignment::from_cigar(const string &a, const string &b, const string &cigar_str)
{
	auto aln = Alignment{ "A", 0, (int)a.size(), "B", 0, (int)b.size(), a, b, "", "", "", {}, {} };
	for (int ci = 0, num = 0; ci < cigar_str.size(); ci++) {
		if (isdigit(cigar_str[ci])) {
			num = 10 * num + (cigar_str[ci] - '0');
		} else if (cigar_str[ci] == ';') {
			continue;
		} else {
			aln.cigar.push_back({cigar_str[ci], num});
			num = 0;
		}
	}
	aln.populate_nice_alignment();
	return aln;
}

string Alignment::cigar_string()
{
	string res;
	for (auto &p: cigar) {
		res += fmt::format("{}{}", p.second, p.first);
	}
	return res;
}

Alignment Alignment::trim() 
{
	auto result = *this;
	while (result.cigar.size()) {
		if (result.cigar[0].first == 'D') {
			result.a = result.a.substr(result.cigar[0].second);
			result.start_a += result.cigar[0].second;
			result.cigar.pop_front();
		} else if (result.cigar[0].first == 'I') {
			result.b = result.b.substr(result.cigar[0].second);
			result.start_b += result.cigar[0].second;
			result.cigar.pop_front();
		} else if (result.cigar.back().first == 'D') {
			result.end_a -= result.cigar.back().second;
			result.a = result.a.substr(0, result.a.size() - result.cigar.back().second);
			result.cigar.pop_back();
		} else if (result.cigar.back().first == 'I') {
			result.end_b -= result.cigar.back().second;
			result.b = result.b.substr(0, result.b.size() - result.cigar.back().second);
			result.cigar.pop_back();
		} else {
			break;
		}
	}

	result.populate_nice_alignment();
	return result;
}

void Alignment::populate_nice_alignment()
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

Alignment::AlignmentError Alignment::calculate_error()
{
	if (!alignment.size()) populate_nice_alignment();

	int gaps = 0, gap_bases = 0, mismatches = 0, matches = 0;
	for (auto &c: cigar) {
		if (c.first != 'M') gaps++, gap_bases += c.second;
	}
	for (int i = 0; i < alignment.size(); i++) {
		if (align_a[i] != '-' && align_b[i] != '-') {
			if (toupper(align_a[i]) == toupper(align_b[i])) {
				matches++;
			} else {
				mismatches++;
			}
		}
	}
	return AlignmentError{gaps, gap_bases, mismatches, matches};
}

string Alignment::print(int width)
{
	if (!alignment.size()) {
		populate_nice_alignment();
	}

	string res;
	int qa = start_a, sa = 0;
	int qb = start_b, sb = 0;

	auto err = calculate_error();
	res += fmt::format(
		"       A: {:>9}..{:<9} (len {:7})    Gaps:       {:5} = {:.0f}% ({})\n"
		"       B: {:>9}..{:<9} (len {:7})    Mismatches: {:5} = {:.0f}%\n"
		"   CIGAR: {}\n",
		start_a, end_a, end_a - start_a, 
		err.gap_bases, err.gap_error(), err.gaps,
		start_b, end_b, end_b - start_b,
		err.mismatches, err.mis_error(),
		cigar_string()
	);
	for (int i = 0; i < alignment.size(); i += width) {
		res += fmt::format(
			"   {:10}: {} {}\n   {:10}  {} {}\n   {:10}: {} {}\n", 
			qa, align_a.substr(i, width),   sa,  
			"", alignment.substr(i, width), i+align_a.substr(i, width).size(),
			qb, align_b.substr(i, width), sb);
		for (auto c: align_a.substr(i, width)) if (c != '-') qa++, sa++;
		for (auto c: align_b.substr(i, width)) if (c != '-') qb++, sb++;
	}
	return res;
}

string Alignment::print_only_alignment(int width)
{
	if (!alignment.size()) {
		populate_nice_alignment();
	}

	if (width == -1) { 
		width = alignment.size();
	}
	string res;
	for (int i = 0; i < alignment.size(); i += width) {
		res += fmt::format(
			"{}\n{}\n{}\n\n", 
			align_a.substr(i, width),  
			alignment.substr(i, width),
			align_b.substr(i, width)
		);
	}
	return res;
}

void Alignment::cigar_from_alignment() 
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

vector<Alignment> Alignment::max_sum(int min_span)
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

	vector<Alignment> results;
	for (auto &h: hits) {
		// translate: a --> b
		results.push_back({
			chr_a, start_a + lookup_a[h.first], start_a + lookup_a[h.second],
			chr_b, start_b + lookup_b[h.first], start_b + lookup_b[h.second],
			a.substr(lookup_a[h.first], lookup_a[h.second] - lookup_a[h.first]),
			b.substr(lookup_b[h.first], lookup_b[h.second] - lookup_b[h.first]),
			align_a.substr(h.first, h.second - h.first),
			align_b.substr(h.first, h.second - h.first),
			alignment.substr(h.first, h.second - h.first),
			{}
		});
		results.back().cigar_from_alignment();
	}

	return results;
}

void Alignment::prepend_cigar(const deque<pair<char, int>> &app)
{
	assert(app.size());
	if (cigar.size() && cigar.front().first == app.back().first) {
		cigar.front().second += app.back().second;
		cigar.insert(cigar.begin(), app.begin(), app.end() - 1);
	} else {
		cigar.insert(cigar.begin(), app.begin(), app.end());
	}
}

void Alignment::append_cigar(const deque<pair<char, int>> &app)
{
	assert(app.size());
	if (cigar.size() && cigar.back().first == app.front().first) {
		cigar.back().second += app.front().second;
		cigar.insert(cigar.end(), next(app.begin()), app.end());
	} else {
		cigar.insert(cigar.end(), app.begin(), app.end());
	}
}

Alignment Alignment::from_anchors(const string &qstr, const string &rstr,
	vector<Hit> &guide,
	const int side) 
{
	// eprn("aligning {} to {}", qstr.size(), rstr.size());

	auto prev = guide.begin();
	Alignment aln = prev->aln;
	for (auto cur = next(prev); cur != guide.end(); cur++) {
		int qs = cur->query_start, qe = cur->query_end;
		int qps = prev->query_start, qpe = prev->query_end;

		int rs = cur->ref_start, re = cur->ref_end;
		int rps = prev->ref_start, rpe = prev->ref_end;

		// eprn("   {}..{}->{}..{}", qps, qpe, rps, rpe);
		// eprn("TO {}..{}->{}..{}", qs, qe, rs, re);

		assert(qpe <= qs);
		assert(rpe <= rs);

		aln.end_a = qe;
		aln.end_b = re;
		aln.a += qstr.substr(qpe, qe - qpe); 
		aln.b += rstr.substr(rpe, re - rpe); 
		
		int qgap = qs - qpe, rgap = rs - rpe;
		if (qgap && rgap) {
			if (qgap <= 1000 && rgap <= 1000) { // "close" hits
				auto gap = align(qstr.substr(qpe, qgap), rstr.substr(rpe, rgap), 5, -4, 40, 1);
				aln.append_cigar(gap.cigar);
			} else { // assume only one part is the gap
				int ma = max(qgap, rgap);
				int mi = min(qgap, rgap);
				auto ma1 = align(qstr.substr(qpe, mi), rstr.substr(rpe, mi), 5, -4, 40, 1);
				ma1.cigar.push_back({qgap == mi ? 'I' : 'D', ma - mi});
				auto ma2 = align(qstr.substr(qs - mi, mi), rstr.substr(rs - mi, mi), 5, -4, 40, 1);
				ma2.cigar.push_front({qgap == mi ? 'I' : 'D', ma - mi});
				aln.append_cigar(ma2.error.error() < ma2.error.error() ? ma2.cigar : ma1.cigar);
			}
		} else if (qgap) {
			aln.append_cigar({{'D', qgap}});	
		} else if (rgap) {
			aln.append_cigar({{'I', rgap}});	
		} 
		// eprn("");
		// assert(qe - qs == re - rs);
		aln.append_cigar(cur->aln.cigar);
		prev = cur;
	}
	// Add end


	int qlo = aln.start_a, qhi = aln.end_a;
	int rlo = aln.start_b, rhi = aln.end_b;
	// eprn("{}", aln.cigar_string());
	// eprn("{}", aln.print());
	// eprn("{} {}\n{} {}", aln.a.size(), aln.a, qhi-qlo, qstr.substr(qlo, qhi - qlo));
	assert(aln.a == qstr.substr(qlo, qhi - qlo));
	assert(aln.b == rstr.substr(rlo, rhi - rlo));
	// eprn("alohaaaa");

	if (side) {
		int qlo_n = max(0, qlo - side);
		int rlo_n = max(0, rlo - side);
		if (qlo - qlo_n && rlo - rlo_n) {
			auto gap = align(
				qstr.substr(qlo_n, qlo - qlo_n), 
				rstr.substr(rlo_n, rlo - rlo_n), 
				5, -4, 40, 1
			);
			// eprn("alignment ok {} vs {}\n{}", qlo-qlo_n, rlo-rlo_n, gap.print());
			gap.trim_front();
			// eprn("trim ok\n{}", gap.print());

			qlo_n = qlo - (gap.end_a - gap.start_a);
			rlo_n = rlo - (gap.end_b - gap.start_b);
			aln.prepend_cigar(gap.cigar);
			aln.a = qstr.substr(qlo_n, qlo - qlo_n) + aln.a;
			aln.b = rstr.substr(rlo_n, rlo - rlo_n) + aln.b;
			aln.start_a = qlo = qlo_n; 
			aln.start_b = rlo = rlo_n;	
		}

		int qhi_n = min(qhi + side, (int)qstr.size());
		int rhi_n = min(rhi + side, (int)rstr.size());
		if (qhi_n - qhi && rhi_n - rhi) {
			auto gap = align(
				qstr.substr(qhi, qhi_n - qhi), 
				rstr.substr(rhi, rhi_n - rhi), 
				5, -4, 40, 1
			);
			gap.trim_back();
			qhi_n = qhi + gap.end_a;
			rhi_n = rhi + gap.end_b;
			aln.append_cigar(gap.cigar);
			aln.a += qstr.substr(qhi, qhi_n - qhi);
			aln.b += rstr.substr(rhi, rhi_n - rhi);
			aln.end_a = qhi = qhi_n;
			aln.end_b = rhi = rhi_n;
		}
	}

	assert(qlo >= 0);
	assert(rlo >= 0);
	assert(qhi <= qstr.size());
	assert(rhi <= rstr.size());
	assert(aln.a == qstr.substr(qlo, qhi - qlo));
	assert(aln.b == rstr.substr(rlo, rhi - rlo));

	// eprn("alohaaaa2");

	aln = aln.trim();
	aln.error = aln.calculate_error();

	// eprn("final {}", aln.cigar_string());
	// eprn("{}", aln.print());
	// cin.get();

	return aln;
}

void Alignment::trim_front() // ABCD -> --CD
{
	int max_score = 0, 
		max_i = a.size();
	int score = 0;
	for (int i = alignment.size() - 1; i >= 0; i--) {
		if (alignment[i] == '|') {
			score += 5;
		} else {
			if (align_a[i] != '-' && align_b[i] != '-') {
				score -= 4;
			} else {
				if (i == alignment.size() - 1 || 
						(align_a[i] == '-' && align_a[i + 1] != '-') || 
						(align_b[i] == '-' && align_b[i + 1] != '-')) 
				{
					score -= 40;
				}
				score -= 1;
			}
		}
		if (score >= max_score) {
			max_score = score, max_i = i;
		}
	}
	if (max_i == a.size())
		return;
	// eprn("max i is {}", max_i);
	// eprn("{}\n{}", cigar_string(), print());
	for (int ci = 0, cur_len = 0; ci < cigar.size(); ci++) {
		// eprn("{} {}", start_a, start_b);
		if (cigar[ci].second + cur_len > max_i) {
			assert(cigar[ci].first == 'M');
			// split this one
			int need = max_i - cur_len;
			cigar[ci].second -= need; 
			for (int cj = 0; cj < ci; cj++)
				cigar.pop_front();
			start_a += need;
			start_b += need;
			break;
		} 
		cur_len += cigar[ci].second;
		if (cigar[ci].first == 'M') {
			start_a += cigar[ci].second; 
			start_b += cigar[ci].second;
		} else if (cigar[ci].first == 'I') {
			start_b += cigar[ci].second;
		} else {
			start_a += cigar[ci].second;
		}
	}
	a = a.substr(start_a, end_a - start_a);
	b = b.substr(start_b, end_b - start_b);
	populate_nice_alignment();
}

void Alignment::trim_back() // ABCD -> AB--
{
	int max_score = 0, 
		max_i = -1;
	int score = 0;
	for (int i = 0; i < alignment.size(); i++) {
		if (alignment[i] == '|') {
			score += 5;
		} else {
			if (align_a[i] != '-' && align_b[i] != '-') {
				score -= 4;
			} else {
				if (i == 0 || 
						(align_a[i] == '-' && align_a[i - 1] != '-') || 
						(align_b[i] == '-' && align_b[i - 1] != '-')) 
				{
					score -= 40;
				}
				score -= 1;
			}
		}
		if (score >= max_score) {
			max_score = score, max_i = i;
		}
	}
	if (max_i == -1)
		return;
	end_a = start_a, end_b = start_b;
	for (int ci = 0, cur_len = 0; ci < cigar.size(); ci++) {
		if (cigar[ci].second + cur_len > max_i) {
			assert(cigar[ci].first == 'M');
			// split this one
			int need = max_i - cur_len;
			cigar[ci].second = need; 
			while (cigar.size() - 1 > ci)
				cigar.pop_back();
			end_a += need;
			end_b += need;
			break;
		} 
		cur_len += cigar[ci].second;
		if (cigar[ci].first == 'M') {
			end_a += cigar[ci].second; 
			end_b += cigar[ci].second;
		} else if (cigar[ci].first == 'I') {
			end_b += cigar[ci].second;
		} else {
			end_a += cigar[ci].second;
		}
	}
	a = a.substr(start_a, end_a - start_a);
	b = b.substr(start_b, end_b - start_b);
	populate_nice_alignment();
}


/******************************************************************************/

Alignment align_helper(const string &qseq, const string &tseq, int sc_mch, int sc_mis, int gapo, int gape, int bandwidth)
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
			bandwidth, -1, // band width; off-diagonal drop-off to stop extension (-1 to disable)
			0, &ez);
		for (int i = 0; i < ez.n_cigar; i++) {
			cigar.push_back({"MDI"[ez.cigar[i] & 0xf], ez.cigar[i] >> 4});
		}
		free(ez.cigar);
	}
	return Alignment{ 
		"A", 0, (int)qseq.size(), 
		"B", 0, (int)tseq.size(), 
		qseq, tseq, "", "", "", cigar 
	};
}

Alignment align(const string &fa, const string &fb, int match, int mismatch, int gap_open, int gap_extend, int bandwidth)
{
	string xa = fa, xb = fb;

	transform(xa.begin(), xa.end(), xa.begin(), align_dna);
	transform(xb.begin(), xb.end(), xb.begin(), align_dna);

	auto a = align_helper(fa, fb, match, mismatch, gap_open, gap_extend, bandwidth);
	a.a = fa, a.b = fb;
	a.populate_nice_alignment();
	a.error = a.calculate_error();
	return a;
}
