/// 786

#pragma once 

#include <map>
#include <vector>
#include <string>

#include "extern/format.h"
#define prn(f, ...)    fmt::print(f "\n", __VA_ARGS__)
#define prnn(...)      fmt::print(__VA_ARGS__)

#define eprn(f, ...)   fmt::print(stderr, f "\n", __VA_ARGS__)
#define eprnn(...)     fmt::print(stderr, __VA_ARGS__)

const int    KMER_SIZE = 14;
// static_assert(KMER_SIZE <= 16, "k-mer space is 32-bit");

const double MAX_GAP_ERROR  = 0.15;
const double MAX_EDIT_ERROR = 0.10;
const double GAP_FREQUENCY  = 0.005;

const int    WINDOW_SIZE    = 16; // <-- Needs to be changed
const int    MIN_READ_SIZE  = 1000;


/// Helper functions

// constexpr auto dna_lookup_init() {
//     using namespace std;
//     array<char, 128> values = {};
//     get<'C'>(values) = 1; get<'c'>(values) = 1;
//     get<'G'>(values) = 2; get<'g'>(values) = 2;
//     get<'T'>(values) = 3; get<'t'>(values) = 3;
//     return values;
// }
// constexpr auto dna_lookup = dna_lookup_init();

static char dna_lookup[128] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,
    3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,
    0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0
};

static char rev_dna[128] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 't', 'N', 'g', 'N', 'N', 'N', 'c', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'a', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};

inline char qdna(char c) 
{
    return dna_lookup[c];
}

inline char rdna(char c) 
{
    return rev_dna[c];
}

template<class X, class Y>
inline bool in_map(const std::map<X, Y> &m, X k)
{
    return m.find(k) != m.end();
}

/// Jaccard helper functions

inline double tau(double error=MAX_EDIT_ERROR)
{
    return (1 - MAX_GAP_ERROR) / (2 * std::exp(KMER_SIZE * error) - 1);
}

inline double j2md(float j)
{
    if (std::fabs(j) < 1e-8) {
        return 100;
    } else if (std::fabs(j - 1) < 1e-8) {
        return 0;
    } else {
        return 100 * (1 - (-1.0 / KMER_SIZE) * std::log(2.0 * j / (1 + j)));
    }
}

int estM(int s, double ci=0.75);
std::vector<std::string> split(const std::string &s, char delim);

double current_time();

