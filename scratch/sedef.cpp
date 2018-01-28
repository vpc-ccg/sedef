/// 786 

#include <string>
using namespace std;

#include "src/common.h"
#include "src/search.h"
#include "src/chain.h"
#include "src/align.h"

class PyHit {
public:
	int qs, qe;
	int rs, re;
	Alignment c;

	// PyHit() = default;
	int query_start() { return qs; }
	int query_end() { return qe; }
	int ref_start() { return rs; }
	int ref_end() { return re; }
	string cigar() { return c.cigar_string(); }
	int alignment_size() { return c.span(); }
	int gaps() { return c.gap_bases();  }
	int mismatches() { return c.mismatches();  }


	bool operator==(const PyHit &q) const {
		return tie(qs, qe, rs, re) == tie(q.qs, q.qe, q.rs, q.re);
	}
};

class PyAligner {
public:
	vector<PyHit> jaccard_align(const string &q, const string &r);
	vector<PyHit> chain_align(const string &q, const string &r);
	vector<PyHit> full_align(const string &q, const string &r);
};

vector<PyHit> PyAligner::jaccard_align(const string &q, const string &r)
{
	shared_ptr<Index> query_hash = make_shared<Index>(make_shared<Sequence>("qry", q));
	shared_ptr<Index> ref_hash = make_shared<Index>(make_shared<Sequence>("ref", r));

	Tree tree;
	vector<PyHit> hits;
	bool iterative = true;
	if (iterative) {
		for (int qi = 0; qi < query_hash->minimizers.size(); qi++) {
			auto &qm = query_hash->minimizers[qi];
			// eprn("search {}", qm.loc);
			if (qm.hash.status != Hash::Status::HAS_UPPERCASE) 
				continue; 
			auto hi = search(qi, query_hash, ref_hash, tree, false);
			for (auto &pp: hi) {
				hits.push_back(PyHit {
					pp.query_start, pp.query_end,
					pp.ref_start, pp.ref_end,
					{}
				});
			}
		}
	} else {
		// Disabled for now
		auto hi = search(0, query_hash, ref_hash, tree, false, max(q.size(), r.size()), false);
		for (auto &pp: hi) {
			hits.push_back({
				pp.query_start, pp.query_end,
				pp.ref_start, pp.ref_end,
				{}
			});
		}
	}
	return hits;
}

vector<PyHit> PyAligner::chain_align(const string &q, const string &r)
{
	extern int DEBUG;
	DEBUG = 0;
	vector<PyHit> hits;
	auto hi = fast_align(q, r);
	for (auto &pp: hi) {
		hits.push_back({
			pp.query_start, pp.query_end,
			pp.ref_start, pp.ref_end,
			pp.aln
		});
	}
	return hits;
}

vector<PyHit> PyAligner::full_align(const string &q, const string &r)
{
	auto aln = Alignment(q, r);
	return vector<PyHit> {{
		0, q.size(),
		0, r.size(),
		aln
	}};
}

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

BOOST_PYTHON_MODULE(sedef)
{
	using namespace boost::python;

    class_<PyHit>("PyHit")
    	.def("cigar", &PyHit::cigar)
    	.def("alignment_size", &PyHit::alignment_size)
    	.def("gaps", &PyHit::gaps)
    	.def("mismatches", &PyHit::mismatches)
    	.def("query_start", &PyHit::query_start)
    	.def("query_end", &PyHit::query_end)
    	.def("ref_start", &PyHit::ref_start)
    	.def("ref_end", &PyHit::ref_end);

	class_<std::vector<PyHit>>("vector")
    	.def(vector_indexing_suite<std::vector<PyHit>, true>());
    
    class_<PyAligner>("PyAligner")
    	.def("jaccard_align", &PyAligner::jaccard_align)
    	.def("chain_align", &PyAligner::chain_align)
    	.def("full_align", &PyAligner::full_align)
    ;
}