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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "align.h"
#include "align_main.h"
#include "common.h"
#include "fasta.h"
#include "chain.h"
#include "hit.h"

using namespace std;

/******************************************************************************/

vector<Hit> merge(vector<Hit> &hits, const int merge_dist) 
{
	vector<Hit> results;

	for (auto &h: hits) {
		assert(h.ref != nullptr);
		assert(h.query != nullptr);
		if (tie(h.query->name, h.query_start, h.query_end) > tie(h.ref->name, h.ref_start, h.ref_end)) {
			swap(h.query->name, h.ref->name);
			swap(h.query_start, h.ref_start);
			swap(h.query_end, h.ref_end);
		}
	}

	sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
		return 
			tie(a.ref->is_rc, a.query->name, a.ref->name, a.query_start, a.ref_start) <
			tie(b.ref->is_rc, b.query->name, b.ref->name, b.query_start, b.ref_start);
	});

	Hit rec, prev;
	int wcount = 0;
	size_t len = 0;
	ssize_t nread;
	multimap<int, Hit> windows;
	for (auto &rec: hits) {
		assert(!rec.query->is_rc);
		if (rec.query->name == rec.ref->name && 
			rec.query_start == rec.ref_start && 
			rec.query_end == rec.ref_end &&
			rec.query->is_rc == rec.ref->is_rc) 
		{
			continue;
		}
		if ((&rec - &hits[0]) == 0) {
			windows.emplace(rec.ref_end, rec);
			prev = rec;
			wcount++;
		} else if (prev.query_end + merge_dist < rec.query_start || 
			prev.query->name != rec.query->name || 
			prev.ref->name != rec.ref->name || 
			prev.ref->is_rc != rec.ref->is_rc) 
		{
			for (auto it: windows) {
				results.push_back(it.second);
			}
			windows.clear();
			windows.emplace(rec.ref_end, rec);
			prev = rec;
			wcount++;
		} else {
			bool needUpdate = 1;
			while (needUpdate) {
				auto start_loc = windows.lower_bound(rec.ref_start - merge_dist);
				needUpdate = 0;
				while (start_loc != windows.end()) {
					if (start_loc->second.query_end + merge_dist < rec.query_start ||
							  start_loc->second.ref_end < rec.ref_start - merge_dist ||
							  start_loc->second.ref_start > rec.ref_end + merge_dist) 
					{
						 start_loc++;
						 continue;
					}
					needUpdate = 1;
					rec.query_end = max(rec.query_end, start_loc->second.query_end);
					rec.ref_end = max(rec.ref_end, start_loc->second.ref_end);
					rec.query_start = min(rec.query_start, start_loc->second.query_start);
					rec.ref_start = min(rec.ref_start, start_loc->second.ref_start);
					windows.erase(start_loc++);
				}
			}
			windows.emplace(rec.ref_end, rec);
		}
		rec.query_end = max(rec.query_end, prev.query_end);
		prev = rec;
	}
	for (auto it: windows)
		results.push_back(it.second);
	return results;
}

/******************************************************************************/

void merge_main(int argc, char **argv)
{
	if (argc < 1) {
		throw fmt::format("Not enough arguments to merge");
	}

	string file = argv[0];

	vector<Hit> hits;
	ifstream fin(file.c_str());
	if (!fin.is_open()) {
		throw fmt::format("BED file {} does not exist", file);
	}

	string s;
	while (getline(fin, s)) {
		Hit h = Hit::from_bed(s);
		hits.push_back(h);
	}
	eprn("Read total {} alignments", hits.size());
	// hits = merge(hits, 0);
	eprn("After merging remaining {} alignments", hits.size());

	sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
	return 
		tie(a.ref->is_rc, a.query->name, a.ref->name, a.query_start, a.ref_start) <
		tie(b.ref->is_rc, b.query->name, b.ref->name, b.query_start, b.ref_start);
	});

	for (auto &h: hits) {
		prn("{}", h.to_bed(false));
	}
}
