/// 786

#include <time.h>
#include <sstream>

#include "common.h"
using namespace std;

vector<string> &split(const string &s, char delim, vector<string> &elems) 
{
    stringstream ss(s);
    string item;
    while(getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim) 
{
    vector<string> elems;
    return split(s, delim, elems);
}

double current_time()
{
    timespec tv;
    clock_gettime(CLOCK_REALTIME, &tv);

    return (tv.tv_sec*1000000000 + tv.tv_nsec) / 1000000000.0;
}