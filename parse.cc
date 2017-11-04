/// 786

#include <bits/stdc++.h>
#include "common.h"
#include "fasta.h"
#include "align.h"

#include "extern/ksw2.h"

using namespace std;

static char rdna[128] = {
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
	'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 't', 'N', 'g', 'N', 'N', 'N', 'c', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
	'N', 'N', 'N', 'N', 'a', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};

// ints are - if rev compl
typedef tuple<string, int, int> loc_t;

string seq(FastaReference &fr, loc_t l) 
{
	int a = get<1>(l);
	int b = get<2>(l);
	bool rc = 0;
	if (a<0) {
		a=-a-2; b=-b-2;
		swap(a,b);
		rc=1;
	}
    string s = fr.getSubSequence(get<0>(l), a, b-a);
    if (rc) {
        reverse(s.begin(), s.end());
        for (auto &c: s) c = rdna[c];
    }
    return s;
}

void prna(loc_t A1, loc_t B1) 
{
	int ll = 0;
	if (get<1>(A1) <= get<1>(B1)) prnn(" A {:10}..<{:^10}> .. {:<10} B ", get<1>(A1), get<1>(B1) - get<1>(A1), get<1>(B1)), ll=get<1>(B1);
	if (get<1>(A1) >  get<1>(B1)) prnn(" B {:10}..<{:^10}> .. {:<10} A ", get<1>(B1), get<1>(A1) - get<1>(B1), get<1>(A1)), ll=get<1>(A1);
	if (get<2>(A1) <= get<2>(B1)) prn(" || <{:^10}> || A {:10} .. <{:^10}> .. {:<10} B ", get<2>(A1) - ll, get<2>(A1), get<2>(B1) - get<2>(A1), get<2>(B1));
	if (get<2>(A1) >  get<2>(B1)) prn(" || <{:^10}> || B {:10} .. <{:^10}> .. {:<10} A ", get<2>(B1) - ll, get<2>(B1), get<2>(A1) - get<2>(B1), get<2>(A1));
}

bool ceq(char a, char b) {
	return toupper(a)==toupper(b);
}

string xalign(string fa, string fb, int MA=5, int MIMA=-4, int GO=40, int GE=-1) ;

bool check_overlap(loc_t x, loc_t y) 
{
	if (get<0>(x) != get<0>(y)) return false;
	if (get<2>(x) < get<1>(y) || get<1>(x) >= get<2>(y)) return false;
	return true;
}

/************************************************************************************************************/
/************************************************************************************************************/
/************************************************************************************************************/
/************************************************************************************************************/
/************************************************************************************************************/

// format: [((i, j) span in edit string, (a, b) start pos in strings)]
vector<pair<pair<int, int>, pair<int, int>>> max_sum(auto al, int THRESHOLD)
{
	const string &aln = get<2>(al);
	const string &a   = get<0>(al);
	const string &b   = get<1>(al);

	vector<pair<pair<int, int>, pair<int, int>>> hits; // start, end

	double score = (aln[0] == '|' ? 1 : -1);
	int max_so_far = score, Mi = 0, Mj = 1;
	int mm = score, mi = 0, mj = 1;

	int AI = 0, BI = 0;
	pair<int, int> mI, MI;

	for (int i = 1; i < aln.size(); i++) {
		score = (aln[i] == '|' ? 1 : -1);
		if (score > mm) {
			if (Mj - Mi >= THRESHOLD && (!hits.size() || hits.back().first != make_pair(Mi, Mj))) {
				hits.push_back({{Mi, Mj}, MI});
			}
			mm = score;
			mi = i, mj = i + 1;
			mI = {AI, BI};
			max_so_far = mm;
			Mi = mi, Mj = mj, MI = mI;
		} else {
			mm += score;
			mj++;
		}
		if (mm > max_so_far) {
			max_so_far = mm;
			Mi = mi, Mj = mj, MI = mI; 
		}
		if (a[i] != '-') AI++;
		if (b[i] != '-') BI++;
	}
	if (Mj - Mi >= THRESHOLD && (!hits.size() || hits.back().first != make_pair(Mi, Mj))) {
		hits.push_back({{Mi, Mj}, MI});
	}
	return hits;
}

void prn_aln (auto result, int qa, int qb, int WD=2000)
{
	for (int i = 0; i < get<2>(result).size(); i += WD) {
		prn("   {:10}: {} {}\n   {:10}  {}\n   {:10}: {}", 
			qa, get<0>(result).substr(i, WD), i+get<0>(result).substr(i, WD).size(),
			"", get<2>(result).substr(i, WD), 
			qb, get<1>(result).substr(i, WD));
		for (auto c: get<0>(result).substr(i, WD)) if (c != '-') qa++;
		for (auto c: get<1>(result).substr(i, WD)) if (c != '-') qb++;
		prnn("\n");
	}
}

pair<double, double> get_err(auto result)
{
	int len = get<2>(result).size();
	int gaps = 0, mismatches = 0;
	for (int i = 0; i < len; i++) {
		if (get<0>(result)[i] == '-' || get<1>(result)[i] == '-')
			gaps++;
		else mismatches += !ceq(get<0>(result)[i], get<1>(result)[i]);
	}
	double gap_err = 100*gaps/double(len);
	double mis_err = 100*mismatches/double(len);
	return {mis_err, gap_err};
}

pair<pair<int, int>, pair<int, int>> trim(string &a, string &b, deque<pair<char, int>> &cigar) 
{
	int fa = 0, fb = 0;
	int ga = 0, gb = 0;
	while (cigar.size()) {
		if (cigar[0].first == 'D') {
			a = a.substr(cigar[0].second);
			prn("  Trim: removing {} at beginning of A1", cigar[0].second);
			fa += cigar[0].second;
			cigar.pop_front();
			continue;
		}
		if (cigar[0].first == 'I') {
			b = b.substr(cigar[0].second);
			prn("  Trim: removing {} at beginning of A2", cigar[0].second);
			fb += cigar[0].second;
			cigar.pop_front();
			continue;
		}

		if (cigar.back().first == 'D') {
			prn("  Trim: removing {} at the back of A1", cigar.back().second);
			ga -= cigar.back().second;
			a = a.substr(0, a.size() - cigar.back().second);
			cigar.pop_back();
			continue;
		}
		if (cigar.back().first == 'I') {
			prn("  Trim: removing {} at the back of A2", cigar.back().second);
			gb -= cigar.back().second;
			b = b.substr(0, b.size() - cigar.back().second);
			cigar.pop_back();
			continue;
		}
		break;
	}
	return {{fa, fb}, {ga, gb}};
}

auto cigar_to_deq(string cigarstr)
{
	deque<pair<char, int>> cigar;
	for (int ci = 0, num = 0; ci < cigarstr.size(); ci++) {
		if (isdigit(cigarstr[ci])) {
			num = 10 * num + (cigarstr[ci] - '0');
		} else if (cigarstr[ci] == ';') {
			continue;
		} else {
			cigar.push_back({cigarstr[ci], num});
			num = 0;
		}
	}
	return cigar;
}

auto check_err(string a, string b, const deque<pair<char, int>> &cigar)
{
	int len = 0;
	int ia = 0, ib = 0;
	tuple<string, string, string> result;
	for (auto &c: cigar) {
		len += c.second;
		for (int i = 0; i < c.second; i++) {
			assert(c.first != 'M' || ia < a.size());
			assert(c.first != 'M' || ib < b.size());
			if (c.first == 'M' && ceq(a[ia], b[ib])) {
				get<2>(result) += "|";
			} else {
				get<2>(result) += " ";
			}
			if (c.first != 'D') get<1>(result) += b[ib++];
			else                get<1>(result) += "-";
			if (c.first != 'I') get<0>(result) += a[ia++];
			else                get<0>(result) += "-";
		}
	}
	auto err = get_err(result);
	return make_pair(err, result);
}

// B1 maps to B2 in B
// In A1, B1 should map to which pos in A2...
// In A2, B2 should map to which pos in A1...
pair<pair<int, int>,pair<int, int>> find_match_in_alignment(auto result, loc_t A1, loc_t A2, loc_t B1, loc_t B2) 
{
	int qa1 = get<1>(A1);
	int qa2 = get<1>(A2); 

	// prn("STA {} ~ {} {}, {} ~ {} {}", qa1, get<1>(B1), get<2>(B1), qa2, get<1>(B2), get<2>(B2));

	pair<int, int> mloc1 = {0, 0};
	pair<int, int> mloc2 = {0, 0};
	for (int i = 0; i < get<2>(result).size(); i++) {
		if (get<0>(result)[i] != '-') qa1++;
		if (get<1>(result)[i] != '-') qa2++;
		
		if (qa1 == get<1>(B1)) mloc1.first  = qa2;
		if (qa1 == get<2>(B1)) mloc1.second  = qa2;

		if (qa2 == get<1>(B2)) mloc2.first  = qa1;
		if (qa2 == get<2>(B2)) mloc2.second  = qa1;
	}
	return make_pair(mloc1, mloc2);
}

bool check (FastaReference &fr, loc_t A1, loc_t A2, loc_t B1, loc_t B2, string cigarstr) 
{
	prnn("Preflight check:\n");

	auto a = seq(fr, A1);
	auto b = seq(fr, A2);
	
	// auto ha = a, hb = b;
	
	auto cigar = cigar_to_deq(cigarstr);
	auto tr = trim(a, b, cigar);
	get<1>(A1) += tr.first.first;  
	get<1>(A2) += tr.first.second; 
	get<2>(A1) -= tr.second.first;
	get<2>(A2) -= tr.second.second;

	prnn("  * 1: "); prna(A1, B1);
	prnn("  * 2: "); prna(A2, B2);		
	
	auto err_aln = check_err(a, b, cigar);

	if (err_aln.first.first + err_aln.first.second <= 25) return true;

	prn("  * len={}; mis={:.1f}; gap={:.1f}; total={:.1f}", 
		get<0>(err_aln.second).size(), 
		err_aln.first.first, err_aln.first.second, 
		err_aln.first.first + err_aln.first.second);

	auto match = find_match_in_alignment(err_aln.second, A1, A2, B1, B2);
	prn("  * >> B1 ({:10} .. {:<10}) and B2 ({:10} .. {:<10}) in A:\n"
		"       B1  {:10} .. {:<10},     B2  {:10} .. {:<10}",
		get<1>(B1), get<2>(B1), get<1>(B2), get<2>(B2),
		match.first.first, match.first.second,
		match.second.first, match.second.second);

	// prn_aln(err_aln.second, get<1>(A1), get<1>(A2), 100);

	// find maximal matches
	auto ncos = max_sum(err_aln.second, 1000); 
	if (ncos.size() == 0) prnn("  !! NO MATCHES !!\n");
	// returns list of max matchings of span at least 750
	// format: [((i, j) span in edit string, (a, b) start pos in strings)]
	// TODO check this for error!!! 	
	for (auto ncop: ncos) {
		auto ea = err_aln;
		auto &nco = ncop.first;
		get<0>(ea.second) = get<0>(ea.second).substr(nco.first, nco.second - nco.first);
		get<1>(ea.second) = get<1>(ea.second).substr(nco.first, nco.second - nco.first);
		get<2>(ea.second) = get<2>(ea.second).substr(nco.first, nco.second - nco.first);
		ea.first = get_err(ea.second);

		if (ea.first.first + ea.first.second > 25) continue;

		auto n_A1 = A1;
		get<1>(n_A1) += ncop.second.first;  // shift A1
		get<2>(n_A1) = get<1>(n_A1);
		for (auto c: get<0>(ea.second)) if (c != '-') get<2>(n_A1)++;

		auto n_A2 = A2;
		get<1>(n_A2) += ncop.second.second; // shift A2
		get<2>(n_A2) = get<1>(n_A2);
		for (auto c: get<1>(ea.second)) if (c != '-') get<2>(n_A2)++;

		prnn("  + 1: "); prna(n_A1, B1);
		prnn("  + 2: "); prna(n_A2, B2);		

		prn("  + len={}; mis={:.1f}; gap={:.1f}; total={:.1f}", 
			get<0>(ea.second).size(), 
			ea.first.first, ea.first.second, 
			ea.first.first + ea.first.second);

		match = find_match_in_alignment(ea.second, n_A1, n_A2, B1, B2);
		prn("  + !! B1 ({:10} .. {:<10}) and B2 ({:10} .. {:<10}) in A:\n"
			"       B1  {:10} .. {:<10},     B2  {:10} .. {:<10}",
			get<1>(B1), get<2>(B1), get<1>(B2), get<2>(B2),
			match.first.first, match.first.second,
			match.second.first, match.second.second);
		prn_aln(ea.second, get<1>(n_A1), get<1>(n_A2));
	}

	// 	//
	// 	//

	// 	auto nerr = get_err(ea.second);
	// 	prn(">=== NEW {}-{} ===<", nco.first, nco.second);
	// 	prn("len={}; mis={:.1f}; gap={:.1f}; total={:.1f}", 
	// 		nco.second - nco.first, nerr.first, nerr.second, nerr.first+nerr.second);
	// 	prn_aln(ea.second, ncos.second[INC].first+get<1>(A1), ncos.second[INC].second+get<1>(A2));

	// }


	// prn_aln(err_aln.second, get<1>(A1), get<1>(A2));


	// prnn("Alt:\n");
	// cigarstr = xalign(a, b, 1,-10,40,1);
	// cigar = cigar_to_deq(cigarstr);
	// err_aln = check_err(a, b, cigar, true, get<1>(A1), get<1>(A2));

	// auto ncos = max_sum(err_aln.second);	

	// auto wa = seq(fr, B1);
	// auto wb = seq(fr, B2);
	// auto wcigstr = xalign(wa, wb);
	// auto wcigar = cigar_to_deq(wcigstr);
	// check_err(wa, wb, wcigar, true, get<1>(B1), get<1>(B2));

	// if (ncos.first.size() == 0) prnn("NO HITS !!!!!!!\n");
	// for (int INC=0;INC<ncos.first.size();INC++) {
	// 	auto &nco = ncos.first[INC];
	// 	auto ea = err_aln;
	// 	get<0>(ea.second) = get<0>(ea.second).substr(nco.first, nco.second - nco.first);
	// 	get<1>(ea.second) = get<1>(ea.second).substr(nco.first, nco.second - nco.first);
	// 	get<2>(ea.second) = get<2>(ea.second).substr(nco.first, nco.second - nco.first);
	// 	auto nerr = get_err(ea.second);
	// 	prn(">=== NEW {}-{} ===<", nco.first, nco.second);
	// 	prn("len={}; mis={:.1f}; gap={:.1f}; total={:.1f}", 
	// 		nco.second - nco.first, nerr.first, nerr.second, nerr.first+nerr.second);
	// 	prn_aln(ea.second, ncos.second[INC].first+get<1>(A1), ncos.second[INC].second+get<1>(A2));
	// }

	return false; // (err.first.first+err.first.second <= 25);
}

loc_t make_loc(string chr, int a, int b, bool rc) 
{
	if (rc) {
		return loc_t(chr, -b, -a);
	} else {
		return loc_t(chr, a, b);
	}
}

void parse_debug(int argc, char **argv)
{
	eprnn("parsiiiing!\n");

	FastaReference fr("data/hg19/hg19.fa");

	ifstream fin(argv[1]);
	string s;
	int _q = 0, i = 0;
	while (getline(fin, s)) {
		prn("{}: {}", ++i, s);
		auto sp = split(s, '\t');

		// Reverse...

		int q = sp[9][0] == '-';
		loc_t B1 = make_loc(sp[0], atoi(sp[1].c_str()), atoi(sp[2].c_str()), false);
		loc_t B2 = make_loc(sp[3], atoi(sp[4].c_str()), atoi(sp[5].c_str()), q);

		q = sp[19][0] == '-';
		loc_t A1 = make_loc(sp[10], atoi(sp[11].c_str()), atoi(sp[12].c_str()), false);
		loc_t A2 = make_loc(sp[13], atoi(sp[14].c_str()), atoi(sp[15].c_str()), q);

		if (!check_overlap(A1, B1)) 
			swap(B1, B2);
		assert(check_overlap(A1, B1));
		assert(check_overlap(A2, B2));

		if (!check(fr, A1, A2, B1, B2, sp[sp.size()-1]))
			_q++;
		if (_q > 10) break;
	}
}

void parse(int argc, char **argv)
{
	FastaReference fr(argv[1]);

	ifstream fin(argv[2]);
	string s;
	while (getline(fin, s)) {
		auto sp = split(s, '\t');

		// Reverse...
		int q = sp[9][0] == '-';
		int q = sp[8][0] == '-';
		loc_t A1 = make_loc(sp[0], atoi(sp[1].c_str()), atoi(sp[2].c_str()), false);
		loc_t A2 = make_loc(sp[3], atoi(sp[4].c_str()), atoi(sp[5].c_str()), q);

		check(fr, A1, A2, sp[sp.size()-1]);
	}
}


// ATAGCTAGCTAGCAT
// --AGCTAcC--GCATA