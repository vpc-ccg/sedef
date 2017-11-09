/// 786

#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <sys/stat.h>
#include <sys/mman.h>

#include "common.h"
#include "fasta.h"
#include "extern/ksw2.h"

using namespace std;

static char align_dna_lookup[128] = {
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
	4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
	3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

// static char align_rev_dna_lookup[128] = {
// 	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
// 	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
// 	4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 2, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
// 	0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 2, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 
// 	4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
// };

string align(const string &tseq, const string &qseq, int sc_mch, int sc_mis, int gapo, int gape)
{
	int8_t a = (int8_t)sc_mch, b = sc_mis < 0 ? (int8_t)sc_mis : (int8_t)(-sc_mis); // a>0 and b<0
	int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
	string cigar;
	const int STEP = 50000;
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
			cigar += fmt::format("{}{}", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
		free(ez.cigar);
		cigar += ";";
	}
	return cigar;
}

string offset_fa(const string &x) 
{
    string s = x;
    for (int i = 0; i < s.size(); i++) {
        s[i] = align_dna_lookup[s[i]];
    }
    return s;
}

string xalign(string fa, string fb, int MA, int MIMA, int GO, int GE) 
{
	for (auto &c: fa) c = align_dna_lookup[c];
	for (auto &c: fb) c = align_dna_lookup[c];
	return align(fa, fb, MA, MIMA, GO, GE);
}

double t_start;
void align(FastaReference &fr, int line, const string &s)
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
        for (auto &c: fb) c = rdna(c);
	}

	string ssa = offset_fa(fa);
	string ssb = offset_fa(fb);
	string aln = align(ssa, ssb, 5, -4, 40, 1);

	eprn("\t{}s\t{}", int(current_time() - t_start), line);
	cout << s << "\t" << aln << "\t" << fa << "\t" << fb << endl;
}

void align(string ref_path, string bed_path, int resume_after) 
{
	int _ = system(">&2 date");

	t_start = current_time();

	FastaReference fr(ref_path);
	ifstream fin(bed_path.c_str());

	int nl = 0;
	string s;
	while (getline(fin, s)) {
		if (nl > resume_after) align(fr, nl, s);
        nl++;
	}

	eprn("Finished bucket {}: read {} lines", string(bed_path), nl);
	eprnn("Done!\n");

	_ = system(">&2 date");
}
