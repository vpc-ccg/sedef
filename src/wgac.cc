/// 786

/******************************************************************************/

#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <thread>
#include <unordered_set>
#include <unordered_map>

#include "common.h"
#include "search.h"
#include "jaccard.h"
#include "fasta.h"
#include "filter.h"
#include "align.h"
#include "aho.h"

using namespace std;

/******************************************************************************/

const int    OFF = 0;
const double MIN_ID = .95;

/******************************************************************************/

double overlap(int sa, int ea, int sb, int eb)
{
	int o = min(ea, eb) - max(sa, sb);
	// eprn("overlap: {} {} / {} {} = {}", sa, ea, sb, eb, o);
	o = max(0, o);

	return o / double(eb - sb);
}

string print_mappings(vector<Hit> &mapping, const string &ca, const string &cb, int la, int lb)
{
	string out;

	sort(mapping.begin(), mapping.end(), [&](Hit pp, Hit pq) {
		auto ta = make_pair(overlap(pp.p, pp.q, OFF, OFF + la), overlap(pp.i, pp.j, OFF, OFF + lb));
		auto tb = make_pair(overlap(pq.p, pq.q, OFF, OFF + la), overlap(pq.i, pq.j, OFF, OFF + lb));
		return ta > tb;
	});
	int cnt = count_if(mapping.begin(), mapping.end(), [](Hit pp){ return pp.reason.substr(0, 2) == "OK"; });
	int prn = 0;
	double prev_score = -1;
	for (auto &pp: mapping) { 
		if (!cnt || pp.reason.substr(0, 2) == "OK") {
			double score = overlap(pp.p, pp.q, OFF, OFF + la) + overlap(pp.i, pp.j, OFF, OFF + lb);
			if (prev_score != -1 && prev_score / score > 2) break;
			out += fmt::format("\n{:.2f} {:.2f} ~ ", 
				overlap(pp.p, pp.q, OFF, OFF + la),
				overlap(pp.i, pp.j, OFF, OFF + lb)) + print_mapping(pp, 0, ca, cb, lb); 
			prev_score = score;
		}
	}
	return out;
}

void check_wgac(string bed_path, string ref_path) 
{
	eprnn("Loading reference... ");
	unordered_map<string, string> ref;
	ifstream ff(ref_path);
	string chr, s;
	bool skip=0;
	while (getline(ff, s)) {
		if (s[0] == '>') {
			chr = s.substr(1);
			// if (chr=="chr1" || chr=="chrX") skip=0;
			// else skip=1;
			if (!skip) eprnn("{} ", chr);
		}
		else if (!skip) ref[chr] += s;
	}
	eprn("done!");

	ifstream fin(bed_path.c_str());
	// getline(fin, s);
	unordered_set<string> seen;
	vector<vector<string>> lines;
	while (getline(fin, s)) {
		auto ss = split(s, '\t');

		if (ss[0][3] == 'U' || ss[0].back() == 'm') continue;
		if (ss[6][3] == 'U' || ss[6].back() == 'm') continue;
		
		// if (ss[0] != "chr1") continue;
		// if (ss[3] != "chr1" && ss[3] != "chrX") continue;

		// if (seen.find(ss[16]) == seen.end()) {
		// 	seen.insert(ss[16]);
			lines.push_back(ss);
		// }
	}
	eprn("IN: {} lines", lines.size());
	// vector<vector<string>> lines = 

	// auto xx =  set<int>{980,990};
	int total = 0, pass = 0, total_fails = 0;
	// #pragma omp parallel for
	for (int si = 0; si < lines.size(); si++) {
		// if (xx.find(si) == xx.end()) continue;
		// if (si != 269) continue;
		eprnn("---- {} ---- ", si);

		auto &ss = lines[si];

		string ca = ss[0]; int sa = atoi(ss[1].c_str()), ea = atoi(ss[2].c_str());
		string cb = ss[3]; int sb = atoi(ss[4].c_str()), eb = atoi(ss[5].c_str());
		bool rb = (ss[9][0] == '_');

		auto refa = ref[ca].substr(sa - OFF, ea - sa);
		auto refb = ref[cb].substr(sb - OFF, eb - sb);
		if (rb) {
			reverse(refb.begin(), refb.end());
			transform(refb.begin(), refb.end(), refb.begin(), rev_dna);
		}

		string out;
		int i_pass = 0, i_total_fails = 0;

		auto aln = alignment_t::from_cigar(refa, refb, ss[11]);
		aln.chr_a = ca; aln.start_a = sa; aln.end_a = ea;
		aln.chr_b = cb; aln.start_b = sb; aln.end_b = eb;
		auto alns = vector<alignment_t>{ aln };
		for (auto &a: alns) {
			auto err = a.calculate_error();
			eprnn(
				"{}\t{}\t{}\t"  
				"{}\t{}\t{}\t"
				"{}\t{:.1f}\t+\t{}\t"
				"{}\t{}\t"
				"SUB:{};{}-{};{}-{}\n",
				aln.chr_a, a.start_a, a.end_a,
				a.chr_b, a.start_b, a.end_b, 
				ss[6], err.error(), ss[9],
				a.alignment.size(), "...",

				a.alignment.size(), // original size
				a.start_a - aln.start_a, a.end_a - aln.start_a, 
				a.start_b - aln.start_b, a.end_b - aln.start_b
			);

			// we need to padd it properly to assure their equality!
			// eprn("{} {}", a.a.size(), a.b.size());
			// eprn("{} {}", refa.size(), refb.size());
			string aa = a.a, ab = a.b;
			if (aa.size() < ab.size()) {
				aa += ref[ca].substr(a.end_a, ab.size() - aa.size());
			}
			if (aa.size() > ab.size()) {
				if (!rb) {
					ab += ref[cb].substr(a.end_b, aa.size() - ab.size());
				} else {
					int real_start = sb + (refb.size() - (a.end_b - a.start_b));
					string rc = ref[cb].substr(real_start - (aa.size() - ab.size()), aa.size() - ab.size());
					reverse(rc.begin(), rc.end());
					transform(rc.begin(), rc.end(), rc.begin(), rev_dna);
					ab += rc;
				}
			}
			assert(aa.size() == ab.size());
			
			Hash ha(aa);
			Hash hb(ab);

			bool success = false;

			TREE_t tree;
			vector<Hit> mappings;

			// 1. Try normal Sedef mapping
			for (int i = OFF; i < ha.seq.size(); i += 250) {
				auto m = search(i, ha, hb, tree, true, MIN_READ_SIZE);
				mappings.insert(mappings.end(), m.begin(), m.end());
				for (auto &pp: m) { 
					if (pp.reason.substr(0, 2) == "OK" 
						&& overlap(pp.p, pp.q, OFF, OFF + aa.size()) >= MIN_ID
						&& overlap(pp.i, pp.j, OFF, OFF + ab.size()) >= MIN_ID)
					{
						success = true;
						break;
					}
				}
				if (success) break;
				// break;
			}
			if (success) {
				eprnn("EXT/OK\n");
				i_pass++;
				continue;
			}

			// eprnn("\n" + a.print());

			eprnn("EXT/FAIL");
			eprnn("{}", print_mappings(mappings, ca, cb, aa.size(), ab.size()));

			// 2. Try full mappings without extension
			mappings = search(OFF, ha, hb, tree, true, max(aa.size(), ab.size()));
			for (auto &pp: mappings) { 
				if (pp.reason.substr(0, 2) == "OK" 
					&& overlap(pp.p, pp.q, OFF, OFF + aa.size()) >= MIN_ID
					&& overlap(pp.i, pp.j, OFF, OFF + ab.size()) >= MIN_ID)
				{
					success = true;
					break;
				}
				if (success) break;
			}
			if (success) {
				eprnn("\nFULL/OK\n");
				continue;
			}
			eprnn("\nFULL/FAIL");
			//eprnn(a.print());
			eprnn("{}", print_mappings(mappings, ca, cb, aa.size(), ab.size()));	
			eprnn("\n");	
			
			i_total_fails++;	
		}


		#pragma omp critical
		{
			eprnn("{}", out);
			total += alns.size();
			total_fails += i_total_fails;
			pass += i_pass;
			// if (i_pass == 0) exit(0);
		}
	}

	eprn("total:     {:>6n}\n"
		  "pass:      {:>6n} ({:.2f})\n"
		  "fail/tot:  {:>6n} ({:.2f})\n"
		  "fail/ext:  {:>6n} ({:.2f})",
		total, 
		pass, pct(pass, total),
		total_fails, pct(total_fails, total),
		total - pass  - total_fails, pct(total - pass  - total_fails, total)
	);
}
