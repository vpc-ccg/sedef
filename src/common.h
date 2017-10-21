/// 786

#include <map>
#include <array>
#include <fmt/format.h>
#include <boost/math/distributions/binomial.hpp>

#define prn(f, ...)    fmt::print(f "\n", __VA_ARGS__)
#define prnn(...)      fmt::print(__VA_ARGS__)

#define eprn(f, ...)    fmt::print(stderr, f "\n", __VA_ARGS__)
#define eprnn(...)      fmt::print(stderr, __VA_ARGS__)

const int    KMER_SIZE = 14;
static_assert(KMER_SIZE <= 16, "k-mer space is 32-bit");

const double ERROR_RATE = 0.10;
const int    WINDOW_SIZE = 16; // <-- Needs to be changed
const int    MIN_READ_SIZE = 1000;

/// Helper functions

constexpr auto dna_lookup_init() {
    using namespace std;
    array<char, 128> values = {};
    get<'C'>(values) = 1; get<'c'>(values) = 1;
    get<'G'>(values) = 2; get<'g'>(values) = 2;
    get<'T'>(values) = 3; get<'t'>(values) = 3;
    return values;
}
constexpr auto dna_lookup = dna_lookup_init();
inline char qdna(char c) 
{
    return dna_lookup[c];
}

template<class X, class Y>
inline bool in_map(const std::map<X, Y> &m, X k)
{
    return m.find(k) != m.end();
}

/// Jaccard helper functions

inline double tau(double error=ERROR_RATE)
{
    return (3.0 / 7) / (2 * std::exp(KMER_SIZE * error) - 1);
}

inline double j2md(float j)
{
    if (std::fabs(j) < 1e-8) {
        return 100;
    } else if (std::fabs(j - 1) < 1e-8) {
        return 0;
    } else {
        return 100 * (1 - (-1.0 / KMER_SIZE) * std::log(2.0 * j/(1 + j)));
    }
}

int estM(int s, double ci=0.75);
