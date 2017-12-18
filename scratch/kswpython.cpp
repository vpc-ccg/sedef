/// 786 

#include <string>
#include <boost/python.hpp>
using namespace std;
#include "extern/ksw2.h"
#include "extern/ksw2_extz2_sse.cc"
#include "extern/format.h"
#include "extern/format.cc"
#include "src/fasta.h"
#include "src/fasta.cc"
#include "src/common.cc"

static char align_dna_lookup[128] = {
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
	4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
	3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

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

void offset_fa(string &x) 
{
    for (auto &c: x) c = align_dna_lookup[c];
}

string kswalign(string sa, string sb, int match=5, int mismatch=-4, int gapopen=40, int gapextend=1)
{
	offset_fa(sa);
	offset_fa(sb);
	string aln = align(sa, sb, match, mismatch, gapopen, gapextend);
	return aln;
}

BOOST_PYTHON_MODULE(kswpython)
{
	using namespace boost::python;
    def("kswalign", kswalign);

    class_<FastaReference>("FastaReference", init<string>())
    	.def("getSubSequence", &FastaReference::getSubSequence)
    ;
}