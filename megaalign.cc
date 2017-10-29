/// 786

#include <bits/stdc++.h>
#include "fmt/fmt/format.h"
#include "fmt/fmt/format.cc"
#include "ksw2/ksw2.h"
#include "ksw2/ksw2_extz2_sse.c"
using namespace std;
#define prn(f, ...)    fmt::print(f "\n", __VA_ARGS__)

static char dna_lookup[128] = {
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,4,
    3,4,4,4,4,4,4,4,4,4,4,4,4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
    4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4
};

string align(string &tseq, string &qseq, int sc_mch, int sc_mis, int gapo, int gape)
{
	for (int i = 0; i < tseq.size(); i++) tseq[i] = dna_lookup[tseq[i]];
	for (int i = 0; i < qseq.size(); i++) qseq[i] = dna_lookup[qseq[i]];
	int a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
	int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
	string cigar;
	const int STEP = 50000;
	for (int SP = 0; SP < min(tseq.size(), qseq.size()); SP += STEP) {
		ksw_extz_t ez;
		ksw_extz2_sse(0, 
			min(STEP, (int)(qseq.size() - SP)), (const uint8_t*)(qseq.c_str() + SP), 
			min(STEP, (int)(tseq.size() - SP)), (const uint8_t*)(tseq.c_str() + SP), 
			5, mat, gapo, gape, -1, -1, 0, &ez);
		
		for (int i = 0; i < ez.n_cigar; ++i) 
			cigar += fmt::format("{}{} ", "MID"[ez.cigar[i]&0xf], ez.cigar[i]>>4);
		free(ez.cigar);
		cigar += ";";
	}
	return cigar;
}

int main(void) {
	std::ios_base::sync_with_stdio(0);
	string na, a, nb, b;
	int i=0;
	while (getline(cin, na)) {
		getline(cin, a);
		getline(cin, nb);
		getline(cin, b);
		fprintf(stderr, "\r %d ~ %d %d          ", ++i, a.size(), b.size());

		// align parameters
		string s = align(a, b, 5, -4, 40, 1);
		for (int i = 0; i < na.size(); i++) 
			if (na[i] == ';') na[i] = '\t';
		prn("{}\t{}", na, s);
	}
	fprintf(stderr, "\r Done!     \n");
	return 0;
}