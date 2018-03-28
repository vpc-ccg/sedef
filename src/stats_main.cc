/// 786 

/******************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <bitset>
#include <unordered_set>
#include <unordered_map>

#include <boost/dynamic_bitset.hpp>

#include "align.h"
#include "common.h"
#include "fasta.h"
#include "hit.h"
#include "merge.h"

using namespace std;

/******************************************************************************/

const int MAX_OK_GAP = 50;
const int MIN_SPLIT_SIZE = 1000;
const int MIN_ALIGNMENT_GAP_SIZE = 100;

/******************************************************************************/

bool subhit(const Hit &hin, int start, int end, Hit &h)
{
	dprn(">subhitting {} {} vs {} {} ", start, end, 0, hin.aln.alignment.size());
	if (end >= hin.aln.alignment.size())
		end = hin.aln.alignment.size();
	if (start >= end) 
		return false;
	h = hin;
	int sa = 0, la = 0, sb = 0, lb = 0;
	for (int i = 0; i < end; i++) {
		if (h.aln.align_a[i] != '-') {
			if (i < start) {
				sa++;
			} else {
				la++;
			}
		}
		if (h.aln.align_b[i] != '-') {
			if (i < start) {
				sb++;
			} else {
				lb++;
			}
		}
	}
	//eprn("{}-{} ({}-{})", hin.query_start, hin.query_end, hin.aln.start_a, hin.aln.end_a);
	//eprn("{}-{} ({}-{})", hin.ref_start, hin.ref_end, hin.aln.start_b, hin.aln.end_b);
	h.aln.align_a = h.aln.align_a.substr(start, end - start);
	h.aln.alignment = h.aln.alignment.substr(start, end - start);
	h.aln.align_b = h.aln.align_b.substr(start, end - start);
	
	h.aln.a = h.aln.a.substr(sa, la);
	h.aln.start_a = 0;
	h.aln.end_a = la;
	
	h.aln.b = h.aln.b.substr(sb, lb);
	h.aln.start_b = 0;
	h.aln.end_b = lb;
	
	h.aln.cigar_from_alignment();
	h.aln.trim_back();
	h.aln.trim_front();
	
	h.query_start += sa;
	h.query_end = h.query_start + la;
	assert(!h.query->is_rc);
	if (h.ref->is_rc) {
		dprn("applying RC");
		h.ref_start = h.ref_end - (lb + sb);
		h.ref_end = h.ref_end - sb;
	} else {
		h.ref_start += sb;
		h.ref_end = h.ref_start + lb;
	}
	return true;
}

vector<Hit> gap_split(Hit h)
{
	struct Gap {
		int start_a, start_b, len_a, len_b;
		int start, len;
	};
	vector<Gap> gaps;
	Gap g {h.aln.start_a, h.aln.start_b, 0, 0, 0, 0};
	for (auto &c: h.aln.cigar) {
		if (c.second && c.first != 'M') {
			if (c.first != 'D') g.len_a = 0, g.len_b = c.second;
			else                g.len_b = 0, g.len_a = c.second;
			g.len = c.second; 
			gaps.push_back(g); 
		}
		if (c.first != 'D') g.start_b += c.second;
		if (c.first != 'I') g.start_a += c.second;
		g.start += c.second;
	}
	sort(gaps.begin(), gaps.end(), [](const Gap &a, const Gap &b) {
		return a.len > b.len;
	});

	Hit hh;
	vector<Hit> hits;
	double gap_score = h.aln.gap_error();
	for (auto &g: gaps) {
		break;
		dprn("--> {} :: a:{} b:{} ... a:{} b:{}", g.len,
			g.start_a - h.aln.start_a, g.start_b - h.aln.start_b,
			h.aln.end_a - (g.start_a + g.len_a), h.aln.end_b - (g.start_b + g.len_b));
		if (g.start_a - h.aln.start_a < MIN_SPLIT_SIZE || g.start_b - h.aln.start_b < MIN_SPLIT_SIZE)
			continue;
		if (h.aln.end_a - (g.start_a + g.len_a) < MIN_SPLIT_SIZE || h.aln.end_b - (g.start_b + g.len_b) < MIN_SPLIT_SIZE)
			continue;

		
		double g_score = pct(g.len, h.aln.error.matches + h.aln.error.gap_bases + h.aln.error.mismatches);

		dprn("{} ~ {}", g_score, g.len);
		if (g_score >= MAX_OK_GAP) {
			// eprn("gap! {} {} {}", g_score, g.len, h.aln.alignment.size());
			dprn(":: Found gap of size {} and score {}\n:: {}\n:: {}\n:: {}", g.len, g_score,
				h.aln.align_a.substr(g.start, g.len),
				h.aln.alignment.substr(g.start, g.len),
				h.aln.align_b.substr(g.start, g.len));

			bool x = subhit(h, 0, g.start, hh);
			assert(x);
			auto r = gap_split(hh);
			for (auto &hx: r) hits.push_back(hx);

			x = subhit(h, g.start + g.len, h.aln.alignment.size(), hh);
			assert(x);
			r = gap_split(hh);
			for (auto &hx: r) hits.push_back(hx);

			return hits;
		}
	}
	if (!hits.size()) hits.push_back(h);
	return hits;
}

void trimlower(Hit &h)
{	
	int sa = 0, sb = 0, i;
	for (i = 0; i < h.aln.align_a.size(); i++) {
		if (isupper(h.aln.align_a[i]) || isupper(h.aln.align_b[i]))
			break;
		sa += bool(h.aln.align_a[i] != '-');
		sb += bool(h.aln.align_b[i] != '-');
	}
	// eprn("sa={} sb={} i={}", sa, sb, i);

	int ea = h.aln.end_a, eb = h.aln.end_b, j;
	for (j = h.aln.align_a.size() - 1; j >= 0; j--) {
		if (isupper(h.aln.align_a[j]) || isupper(h.aln.align_b[j]))
			break;
		ea -= bool(h.aln.align_a[j] != '-');
		eb -= bool(h.aln.align_b[j] != '-');
	}
	j++;
	// eprn("sa={} sb={} i={}", ea, eb, j);

	h.aln.align_a = h.aln.align_a.substr(i, j - i);
	h.aln.alignment = h.aln.alignment.substr(i, j - i);
	h.aln.align_b = h.aln.align_b.substr(i, j - i);
	
	h.aln.a = h.aln.a.substr(sa, ea - sa);
	h.aln.start_a = sa;
	h.aln.end_a = ea;
	
	h.aln.b = h.aln.b.substr(sb, eb - sb);
	h.aln.start_b = sb;
	h.aln.end_b = eb;
	
	h.aln.cigar_from_alignment();
	// h.aln.trim_back();
	// h.aln.trim_front();
	
	h.query_start += sa;
	h.query_end = h.query_start + (ea - sa);
	assert(!h.query->is_rc);
	if (h.ref->is_rc) {
		dprn("applying RC");
		h.ref_start = h.ref_end - eb;
		h.ref_end = h.ref_end - sb;
	} else {
		h.ref_start += sb;
		h.ref_end = h.ref_start + (eb - sb);
	}
}

vector<Hit> split_alignment(Hit h) 
{
	// Find all alignments!!!
	vector<Hit> hits;

	// Find stretch of Ns
	int prev_an = 0, prev_bn = 0;
	int hit_begin = 0;
	Hit hh;
	for (int i = 0; i < h.aln.alignment.size(); i++) {
		if (toupper(h.aln.align_a[i]) == 'N') {
			prev_an++;
		} else {
			if (prev_an >= MIN_ALIGNMENT_GAP_SIZE) {
				dprn(":: Found assembly gap of size {}\n:: {}\n:: {}\n:: {}", prev_an,
					h.aln.align_a.substr(i - prev_an, prev_an),
					h.aln.alignment.substr(i - prev_an, prev_an),
					h.aln.align_b.substr(i - prev_an, prev_an));
				if (subhit(h, hit_begin, i - prev_an, hh))
					hits.push_back(hh);
				hit_begin = i;
			}
			prev_an = 0;
		}
		if (toupper(h.aln.align_b[i]) == 'N') {
			prev_bn++; 
		} else {
			if (prev_bn >= MIN_ALIGNMENT_GAP_SIZE) {
				dprn(":: Found assembly gap of size {}\n:: {}\n:: {}\n:: {}", prev_bn,
					h.aln.align_a.substr(i - prev_bn, prev_bn),
					h.aln.alignment.substr(i - prev_bn, prev_bn),
					h.aln.align_b.substr(i - prev_bn, prev_bn));
				if (subhit(h, hit_begin, i - prev_bn, hh)) 
					hits.push_back(hh);
				hit_begin = i;
			}
			prev_bn = 0;
		}
	}
	if (!hit_begin)
		hits.push_back(h);
	else if (subhit(h, hit_begin, h.aln.alignment.size(), hh))
		hits.push_back(hh);

	// Find gaps
	vector<Hit> hits_final;
	for (auto &h: hits) {
		auto hh = gap_split(h);
		// dprn("first round: {} ret", hh.size());
		for (auto &hhh: hh) {
			hits_final.push_back(hhh);
		}
	}
	return hits_final;
}

void process(Hit hs, string cigar, FastaReference &fr)
{
	// auto ita = reference.find(hs.query->name);
	// assert(ita != reference.end());
	// auto itb = reference.find(hs.ref->name);
	// assert(itb != reference.end());

	string fa = fr.get_sequence(hs.query->name, hs.query_start, &hs.query_end);
	string fb = fr.get_sequence(hs.ref->name, hs.ref_start, &hs.ref_end);
	// if (hs.query_end > ita->second.size()) 
	// 	hs.query_end = ita->second.size();
	// string fa = ita->second.substr(hs.query_start, hs.query_end - hs.query_start);
	// if (hs.ref_end > itb->second.size()) 
	// 	hs.ref_end = itb->second.size();
	// string fb = itb->second.substr(hs.ref_start, hs.ref_end - hs.ref_start);
	assert(!hs.query->is_rc); 
	if (hs.query->is_rc) {
		fa = rc(fa);
	}
	if (hs.ref->is_rc) {
		fb = rc(fb);
	}
	assert(cigar.size());
	hs.aln = Alignment(fa, fb, cigar);

	auto hs_split = split_alignment(hs);
	for (auto &h: hs_split) {
		// eprn("{}", h.aln.print(80));
		// trimlower(h);
		if (h.aln.alignment.size() < 900)
			continue;
		// if (!h.aln.alignment.size())
			// continue;
		// eprn("alnd");
		// eprn("{}", h.aln.print(80));

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
		int uppercaseMatches = 0;

		for (int i = 0; i < align_length; i++) {
			char a = toupper(h.aln.align_a[i]);
			char b = toupper(h.aln.align_b[i]);
			indel_a += a == '-';
			indel_b += b == '-';
			matchB += a != '-' && a == b;
			uppercaseA += (h.aln.align_a[i] != '-' && toupper(h.aln.align_a[i]) != 'N' 
				&& isupper(h.aln.align_a[i]));
			uppercaseB += (h.aln.align_b[i] != '-' && toupper(h.aln.align_b[i]) != 'N' 
				&& isupper(h.aln.align_b[i]));
			if (a != '-' && b != '-') {
				alignB += 1;
				if (a != b) {
					mismatchB += 1;
					if (a == 'A' || a == 'G') {
						transitionsB += b == 'A' || b == 'G';
						transversionsB += !(b == 'A' || b == 'G');
					} else {
						transitionsB += b == 'C' || b == 'T';
						transversionsB += !(b == 'C' || b == 'T');
					}
				} else if (isupper(h.aln.align_a[i]) && isupper(h.aln.align_b[i])) {
					uppercaseMatches++;
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

		// Split large gaps?
		if (uppercaseA >= 100 
			&& uppercaseB >= 100 
			&& !too_big_overlap 
			&& errorScaled <= .50
		 	&& uppercaseMatches >= 100
		 	)
		{
			string l = h.to_bed(false);
			
			h.name = "S";
			h.comment = ""; 
			#pragma omp critical
			prn("{}\t"
				"{}\t{}\t{}\t{}\t{}\t{}\t"
				"{}\t{}\t{}\t{}\t{}\t{}\t"
				"{}\t{}\t{}\t"
				"{}\t{}\t{}\t{}\t"
				"{}\t{}", 
				h.to_bed(false, false), // 1-13
				align_length, // 14
				indel_a, indel_b, // 15-16
				alignB, matchB, mismatchB, // 17-19
				transitionsB, transversionsB, // 20-21
				fracMatch, fracMatchIndel, // 22-23
				jcK, k2K, // 24-25
				h.aln.gaps(), // 26
				uppercaseA, uppercaseB, uppercaseMatches, // 27-29
				h.aln.matches(), h.aln.mismatches(), h.aln.gaps(), h.aln.gap_bases(), // 30-33
				h.aln.cigar_string() // 34
			);
		}
	}
}

void stats(const string &ref_path, const string &bed_path) 
{
	FastaReference fr(ref_path);
	ifstream fin(bed_path.c_str());
	if (!fin.is_open()) {
		throw fmt::format("BED file {} does not exist", bed_path);
	}

	// string sq = "chr7	14840644	14858402	chr7	16454676	16460667	S248049	72.9	+	+	17758	18059	108M1I210M10I76M1I307M1I568M1D45M1D475M1I1677D99M1I106M1D480M2I88M28D61M7D46M7D156M1I41M4D24M17D96M32I43M9D51M4I75M3D21M23I13M13I1M172I57M7D41M1D184M1I23M10D78M1D125M1D133M5D27M1I21M1D17M1I119M7D64M1D139M1I30M1D13M2I61M1D218M1I30M29I10277D352M1I587M1I165M1I16M";
	// // string sq = "chr17	13718022	13742917	chr17	20899063	20929192	S382400	53.2	+	-	30129	36658	308M1D114M6I179M4D208M1I71M4D9M1I288M2I409M1D103M1D200M3D31M16I149M6I544M1I513M2I254M6D99M2I21M2D2731I166M47D100M1I286M4D51M2I31M1D286M10D42M5I111M7I258M1I125M1I99M4D78M22D136M22I50M1278I97M1I273M3D138M29D368M2I28M1I170M2I1072M3D187M845I158M11I174M1D235M1I488M794D312M25D74M25I6704I42M1D104M1I19M1D686M1I84M1I475M3I158M5135D75M1I241M1D155M5D253M111D26M58D15M151D31M63D27M12I31M4D78M2I75M3I21M1D568M1D1618M1D87M4D988M1I56M1D330M42I154M1D22M2I60M1D22M1I438M1I129M1I172M6I46M1I164M17D492M1D175M1D22M3D25M1D145M1I24M3I322M1D131M2I363M1I124M";
	// string cigar;
	// Hit hs = Hit::from_bed(sq, &cigar);
	// string fa, fb;
	// fa = fr.get_sequence(hs.query->name, hs.query_start, &hs.query_end);
	// fb = fr.get_sequence(hs.ref->name, hs.ref_start, &hs.ref_end);
	// assert(!hs.query->is_rc); 
	// if (hs.query->is_rc) {
	// 	fa = rc(fa);
	// }
	// if (hs.ref->is_rc) {
	// 	fb = rc(fb);
	// }
	// assert(cigar.size());
	// hs.aln = Alignment(fa, fb, cigar);
	// dprn(">> IN: {}\t{}", hs.to_bed(0,0), hs.aln.cigar_string());

	// auto hh = split_alignment(hs);
	// for (auto &h: hh) {
	// 	dprn(">> OU: {}\t{}", h.to_bed(0,0), h.aln.cigar_string());
	// 	dprn("{}", h.aln.print(80));
	// }
	// exit(0);

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
			for (auto &c: cigar) {
				if (c == 'I') c = 'D';
				else if (c == 'D') c = 'I';
			}
			// h.aln.swap();
		}
		hits.push_back({h, cigar});
	}
	sort(hits.begin(), hits.end(), [](const pair<Hit, string> &a, const pair<Hit, string> &b) {
		return 
			tie(a.first.ref->is_rc, a.first.query->name, a.first.ref->name, a.first.query_start, a.first.ref_start) <
			tie(b.first.ref->is_rc, b.first.query->name, b.first.ref->name, b.first.query_start, b.first.ref_start);
	});

	eprn("\nRead {:n} hits", hits.size());
	int hit_count = 0, out_count = 0;
	string prev;

	#pragma omp parallel for
	for (auto hsi = 0; hsi < hits.size(); hsi++) {
		process(hits[hsi].first, hits[hsi].second, fr);
		#pragma omp critical
		eprnn("\rProcessed hit {:n}", ++hit_count);
	}
	eprn("\nRead {:n} hits, wrote {:n} SDs", hit_count, out_count);
	eprn("\nDone!");
}

/******************************************************************************/

void get_differences(const string &ref_path, const string &bed_path,
	const string &wgac_path)
{
	map<string, boost::dynamic_bitset<>> sedef;
	map<string, boost::dynamic_bitset<>> wgac;

	string s;
	ifstream fin(bed_path);
	while (getline(fin, s)) {
		string cigar;
		Hit h = Hit::from_bed(s, &cigar);

		auto c1 = fmt::format("{}", h.query->name, "+-"[h.query->is_rc]);
		auto c2 = fmt::format("{}", h.ref->name, "+-"[h.ref->is_rc]);
		if (sedef.find(c1) == sedef.end()) {
			sedef[c1] = boost::dynamic_bitset<>(250 * MB);
		}
		if (sedef.find(c2) == sedef.end()) {
			sedef[c2] = boost::dynamic_bitset<>(250 * MB);
		}
		for (int i = h.query_start; i < h.query_end; i++) 
			sedef[c1].set(i);
		for (int i = h.ref_start; i < h.ref_end; i++) 
			sedef[c2].set(i);
	}

	eprn("SEDEF reading done!");

	ifstream fiw(wgac_path);
	getline(fiw, s);
	unordered_set<string> seen;
	while (getline(fiw, s)) {
		Hit h = Hit::from_wgac(s);
		auto c1 = fmt::format("{}", h.query->name, "+-"[h.query->is_rc]);
		auto c2 = fmt::format("{}", h.ref->name, "+-"[h.ref->is_rc]);		
		if (c1.size() > 6 || c2.size() > 6)
			continue;
		
		if (seen.find(h.name) == seen.end()) {
			seen.insert(h.name);
			if (wgac.find(c1) == wgac.end()) {
				wgac[c1] = boost::dynamic_bitset<>(250 * MB);
			}
			if (wgac.find(c2) == wgac.end()) {
				wgac[c2] = boost::dynamic_bitset<>(250 * MB);
			}
			for (int i = h.query_start; i < h.query_end; i++) 
				wgac[c1].set(i);
			for (int i = h.ref_start; i < h.ref_end; i++) 
				wgac[c2].set(i);
		}
	}

	eprn("WGAC reading done!");

	FastaReference fr(ref_path);

	int intersect = 0, wgac_only = 0, wgac_span = 0, sedef_only = 0, sedef_span = 0;

	int sedef_extra_upper = 0;
	int miss_upper = 0;

	for (auto &p: sedef) {
		auto &s = p.second;
		auto &w = wgac[p.first];

		auto seq = fr.get_sequence(p.first);

		for (int i = 0; i < seq.size(); i++) {
			if ((s[i] & (~w[i])) && isupper(seq[i]) && toupper(seq[i]) != 'N') {
				sedef_extra_upper++;
			}
			if ((w[i] & (~s[i])) && isupper(seq[i]) && toupper(seq[i]) != 'N') {
				miss_upper++;
			}
		}

		intersect += (s & w).count();
		wgac_only += (w & (~s)).count();
		sedef_only += (s & (~w)).count();
		sedef_span += s.count();
		wgac_span += w.count();
	}
	
	eprn("SEDEF: spans              {:12n}\n"
		 "       unique             {:12n}\n"
		 "       unique (uppercase) {:12n}\n"
		 "       misses             {:12n}\n"
		 "       misses (uppercase) {:12n}\n"
		 "WGAC:  spans              {:12n}\n"
		 "       intersects         {:12n}", 
		 sedef_span, sedef_only, sedef_extra_upper, wgac_only, 
		 miss_upper, wgac_span, intersect);
}

/******************************************************************************/

void stats_main(int argc, char **argv)
{
	if (argc < 3) {
		throw fmt::format("Not enough arguments to stats");
	}

	string command = argv[0];
	if (command == "generate") {
		stats(argv[1], argv[2]);
	} else if (command == "diff") {
		get_differences(argv[1], argv[2], argv[3]);
	} else {
		throw fmt::format("Unknown stats command");
	}
}
