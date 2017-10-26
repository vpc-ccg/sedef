/// 786

#include "common.h"
#include <boost/math/distributions/binomial.hpp>

int estM(int s, double ci)
{
    using namespace boost::math;

    int i = ceil(s * ::tau());
    for (; i >= 0; i--)
    {
        int x = quantile(complement(binomial(s, double(i) / s), (1.0 - ci) / 2));
        if (j2md(double(x) / s) >= MAX_GAP_ERROR)
            break;
    }
    return i ? i - 1 : 0;
}


