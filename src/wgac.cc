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
#include "fasta.h"
#include "filter.h"
#include "align.h"
#include "aho.h"
#include "sliding.h"

using namespace std;

/******************************************************************************/

const double MIN_ID = .95; // Minimum overlap percentage to call successful WGAC hit

/******************************************************************************/

double overlap(int sa, int ea, int sb, int eb)
{
	int o = max(0, min(ea, eb) - max(sa, sb));
	return o / double(eb - sb);
}

/******************************************************************************/

vector<Hit> read_wgac(string ref_path, string tab_path, bool is_wgac, string chr = "")
{
	eprnn("Loading reference... ");
	FastaReference fr(ref_path);
	map<pair<string, bool>, shared_ptr<Sequence>> ref;

	string s;
	ifstream fin(tab_path.c_str());
	if (!fin.is_open()) {
		throw fmt::format("WGAC TAB file {} does not exist", tab_path);
	}
	if (is_wgac) {
		getline(fin, s); // header
	}

	unordered_set<string> seen;
	vector<Hit> hits;
	while (getline(fin, s)) {
		Hit hit;
		if (is_wgac) {
			hit = Hit::from_wgac(s);		
		} else {
			hit = Hit::from_bed(s);		
		}

		string chrq = hit.query->name, 
		       chrr = hit.ref->name;
		assert(!hit.query->is_rc);
		bool rc = hit.ref->is_rc;

		if (chrq.size() > 5 || chrr.size() > 5)
			continue;
		if (chr != "" && (chrq != chr || chrr != chr))
			continue;

		if (ref.find({chrq, false}) == ref.end()) 
			ref[{chrq, false}] = make_shared<Sequence>(chrq, fr.get_sequence(chrq));
		if (ref.find({chrr, rc}) == ref.end()) 
			ref[{chrr, rc}] = make_shared<Sequence>(chrr, fr.get_sequence(chrr), rc);
		hit.query = ref[{chrq, false}];
		hit.ref = ref[{chrr, rc}];

		if (hit.ref->is_rc) {
			swap(hit.ref_start, hit.ref_end);
			hit.ref_start = hit.ref->seq.size() - hit.ref_start;
			hit.ref_end = hit.ref->seq.size() - hit.ref_end;
		}
		if (!is_wgac || seen.find(hit.name) == seen.end()) {
			seen.insert(hit.name);
			hits.push_back(hit);
		}
	}
	eprn("Loaded {} hits", hits.size());

	return hits;
}

/******************************************************************************/

//38781571	38786875
// 12517691
void align_wgac(string ref_path, string tab_path)
{
	auto hits = read_wgac(ref_path, tab_path, /*is_wgac*/ true, "chr22");
	vector<Hit> fast_align(const string &sa, const string &sb);
		
	// #pragma omp parallel for
	for (int si = 0; si < hits.size(); si++) {
		auto &hit = hits[si];
		// if (hit.name != "align_both/0015/both076775") continue;
		
		eprn("{}\n{}", string(100, '*'), hit.to_bed());

		auto refq = hit.query->seq.substr(hit.query_start, hit.query_end - hit.query_start);
		auto refr = hit.ref->seq.substr(hit.ref_start, hit.ref_end - hit.ref_start);
		auto hits = fast_align(refq, refr);
		assert(hits.size() != 0);
		eprn("woohoo, size={}", hits.size());
		// eprn("{}", refq);
		// eprn("{}", refr);
		
		// #pragma omp critical 
		auto best = hits.front();
		best.comment = "http://humanparalogy.gs.washington.edu/build37/" + hit.name + " ";

		double offq = best.query_start;
		offq += refq.size() - best.query_end;
		
		double offr = best.ref_start;
		offr += refr.size() - best.ref_end;
		
		eprn("{:3.0f} {:3.0f} ~> {:3.0f} {:3.0f} |> {} |> {}", offq*100/refq.size(), offr*100/refr.size(), offq, offr, best.to_bed(), hit.comment);

		// if (best.ref->is_rc && si > 10) {
		// 	// prn("{}", best.aln.print(80));
		// 	// break;
		// }
	}
}

void check_wgac(string ref_path, string bed_path) 
{
	#define parprnn(...) out += fmt::format(__VA_ARGS__)
	auto print_mappings = [](vector<Hit> &mapping, int la, int lb, string end = "\n") {
		string out;
		sort(mapping.begin(), mapping.end(), [&](const Hit &pp, const Hit &pq) {
			auto ta = make_pair(overlap(pp.query_start, pp.query_end, 0, la), overlap(pp.ref_start, pp.ref_end, 0, lb));
			auto tb = make_pair(overlap(pq.query_start, pq.query_end, 0, la), overlap(pq.ref_start, pq.ref_end, 0, lb));
			return ta > tb;
		});
		int cnt = count_if(mapping.begin(), mapping.end(), [](Hit pp){ return pp.comment.substr(0, 2) == "OK"; });
		int prn = 0;
		double prev_score = -1;
		for (auto &pp: mapping) { 
			if (!cnt || pp.comment.substr(0, 2) == "OK") {
				double score = overlap(pp.query_start, pp.query_end, 0, la) 
					+ overlap(pp.ref_start, pp.ref_end, 0, lb);
				if (prev_score != -1 && prev_score / score > 2) break;
				out += fmt::format("   {:.2f} {:.2f} # ", 
					overlap(pp.query_start, pp.query_end, 0, la),
					overlap(pp.ref_start, pp.ref_end, 0, lb)) + pp.to_bed() + end; 
				prev_score = score;
			}
		}
		return out;
	};

	auto hits = read_wgac(ref_path, bed_path, false);
	int total = 0, pass = 0, total_fails = 0;

	// #pragma omp parallel for
	for (int si = 0; si < 10 /*hits.size()*/; si++) {
		auto &hit = hits[si];
		string out;

		auto refq = hit.query->seq.substr(hit.query_start, hit.query_end - hit.query_start);
		auto refr = hit.ref->seq.substr(hit.ref_start, hit.ref_end - hit.ref_start);
		throw "please don't";
		//hit.aln = Alignment::from_cigar(refq, refr, CIGAR__);
		
		
		parprnn("{}\t", hit.to_bed());

		auto query_hash = make_shared<Index>(make_shared<Sequence>("A", refq)), 
		     ref_hash = make_shared<Index>(make_shared<Sequence>("B", refr));
		
		vector<Hit> mappings;
		Tree tree;
		bool success = false;

	// 1. Try normal Sedef mapping
		int i_pass = 0, i_total_fails = 0;
		for (int i = 0; i < ref_hash->minimizers.size(); i++) {
			auto &qm = ref_hash->minimizers[i];
			if (qm.hash.status != Hash::Status::HAS_UPPERCASE) 
				continue;
			auto m = search(i, query_hash, ref_hash, tree, /*same_genome=*/ false);
			mappings.insert(mappings.end(), m.begin(), m.end());
			for (auto &pp: m) { 
				if (pp.comment.substr(0, 2) == "OK" && 
					overlap(pp.query_start, pp.query_end, 0, refq.size()) >= MIN_ID && 
					overlap(pp.ref_start, pp.ref_end, 0, refr.size()) >= MIN_ID)
				{
					success = true;
					break;
				}
			}
			if (success) {
				parprnn("EXT/OK\n");
				i_pass++;
				break;
			}
		}
		if (!success) {
	// 2. Try full mappings without extension
			parprnn("\n  EXT/FAIL\n");
			parprnn("{}", print_mappings(mappings, refq.size(), refr.size()));

			mappings = search(0, query_hash, ref_hash, tree, 
				/*same_genome=*/ false, 
				/*init_len=*/ max(refr.size(), refq.size()), 
				/*allow_extend=*/ false, 
				/*report_fails=*/ true);
			for (auto &pp: mappings) { 
				if (pp.comment.substr(0, 2) == "OK" 
					&& overlap(pp.query_start, pp.query_end, 0, refq.size()) >= MIN_ID
					&& overlap(pp.ref_start, pp.ref_end, 0, refr.size()) >= MIN_ID)
				{
					success = true;
					break;
				}
				if (success) {
					parprnn("  FULL/OK\n");
					break;
				}	
			}
			if (!success) {
				parprnn("  FULL/FAIL\n");
				// parprnn(a.print(a.alignment.size()));
				parprnn("{}", print_mappings(mappings, refq.size(), refr.size()));
				i_total_fails++;	
			}
		}

		#pragma omp critical
		{
			if (i_pass == 0) prnn("{}", out);
			total++;
			total_fails += i_total_fails;
			pass += i_pass;
		}
	}

	eprn("Total:        {:>6n}\n"
	     "Pass:         {:>6n} ({:.2f})\n"
	     "Fail [extn]:  {:>6n} ({:.2f})"
	     "Fail [full]:  {:>6n} ({:.2f})\n",
		total, 
		pass, pct(pass, total),
		total_fails, pct(total_fails, total),
		total - pass  - total_fails, pct(total - pass  - total_fails, total)
	);
}

/******************************************************************************/

void wgac_main (int argc, char **argv)
{
	// aho = make_shared<AHOAutomata>();

	if (argc < 3) {
		throw fmt::format("Not enough arguments to wgac");
	}

	string command = argv[0];
	if (command == "align") {
		align_wgac(argv[1], argv[2]);
	} else if (command == "check") {
		check_wgac(argv[1], argv[2]);
	} else {
		throw fmt::format("Unknown wgac command");
	}
}

