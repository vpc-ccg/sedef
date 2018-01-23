/// 786

/******************************************************************************/

#include <time.h>
#include <sstream>
#include <unordered_map>

#include <boost/math/tools/roots.hpp>
#include <boost/math/distributions/binomial.hpp>

#include "common.h"
using namespace std;

/******************************************************************************/

vector<string> split(const string &s, char delim) 
{
	vector<string> elems;
	stringstream ss(s);
	string item;
	while(getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

string rc(const string &s)
{
	auto r = s;
	reverse(r.begin(), r.end());
	transform(r.begin(), r.end(), r.begin(), rev_dna);
	return r;
}

/******************************************************************************/

double tau(double edit_error)
{
    double gap_error = std::min(1.0, ERROR_RATIO * edit_error);
    double a = (1 - gap_error) / (1 + gap_error);
    double b = 1 / (2 * std::exp(KMER_SIZE * edit_error) - 1);
    return a * b;
}

double solve_inverse_jaccard(int j)
{
	if (j == 0) return 1;
	if (j == 1) return 0;
	return boost::math::tools::newton_raphson_iterate([j](double d){
		double E = exp(d * KMER_SIZE);
		return make_tuple(
			((1 - d * ERROR_RATIO) / (1 + d * ERROR_RATIO)) * (1.0 / (2 * E - 1)) - j,
			2 * (- KMER_SIZE * E + ERROR_RATIO - 2 * ERROR_RATIO * E + E * KMER_SIZE * pow(d * ERROR_RATIO, 2)) /
				pow((2 * E - 1) * (1 + d * ERROR_RATIO), 2)
		);
	}, 0.10, 0.0, 1.0, numeric_limits<double>::digits);
}

int relaxed_jaccard_estimate(int s)
{
	static unordered_map<int, int> mm;

	double result = -1;
	// #pragma omp critical
	{
		auto it = mm.find(s);
		if (it != mm.end()) result = it->second;
	}
	if (result != -1) return result;

	using namespace boost::math;
	const double CI = 0.75;
	const double Q2 = (1.0 - CI) / 2; // one side interval probability

	result = ceil(s * tau());
	for (; result >= 0; result--) {        
		double d = solve_inverse_jaccard(result / s); // returns edit error
		// eprn("[in={}] inverse for {} is {} ~ {}", s, result/s, d, tau(d));
		double x = quantile(complement(binomial(s, tau(d)), Q2)); // inverse binomial 
		double low_d = solve_inverse_jaccard(x / s);
		if (100 * (1 - low_d) < MAX_EDIT_ERROR) {
			result++; 
			break;
		}
	}
	result = max(result, 0.0);
	// #pragma omp critical
	{
		mm[s] = result;
	}
	return result;
}
