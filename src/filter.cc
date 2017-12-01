/// 786
/// Based on http://www.biorxiv.org/content/biorxiv/early/2017/03/24/103812.full.pdf

#include <list>
#include <vector>
#include <string>
#include <cmath>
#include "common.h"
#include "filter.h"
using namespace std;

int64_t QGRAM_NORMAL_FAILED = 0;
int64_t OTHER_FAILED = 0;
int64_t CORE_FAILED = 0;

AHOAutomata *aho = NULL;
vector<int> qgram_p, qgram_r;

const int MIN_UPPERCASE = 500;

inline int min_qgram(int l, int q) 
{
    return l * (1 - MAX_GAP_ERROR - q * MAX_EDIT_ERROR) - (GAP_FREQUENCY * l + 1) * (q - 1);
}

bool filter(const string &q, int q_pos, int q_len, const string &r, int r_pos, int r_len) 
{
	if (!qgram_p.size()) {
		qgram_p = vector<int>(1 << (2 * 8), 0);
		qgram_r = vector<int>(1 << (2 * 8), 0);
		aho = new AHOAutomata();
	}

    int q_up = 0; for (int i = 0; i < q_len; i++) q_up += isupper(q[q_pos + i]); 
    int r_up = 0; for (int i = 0; i < r_len; i++) r_up += isupper(r[r_pos + i]);

    if (q_up < MIN_UPPERCASE || r_up < MIN_UPPERCASE) {
        OTHER_FAILED++;
        return false;
    }
    
    int maxlen = max(q_len, r_len);
    int QG = 5;
    uint64_t QSZ = (1 << (2 * QG)); 
    uint64_t MASK = QSZ - 1;

    int minqg = min_qgram(maxlen, QG);
    assert(minqg >= 10);

    assert(q_pos + q_len < q.size());
    assert(r_pos + r_len < r.size());
    
    for (uint64_t qi = q_pos, qgram = 0; qi < q_pos + q_len; qi++) {
        qgram = ((qgram << 2) | qdna(q[qi])) & MASK;
        if (qi - q_pos >= QG - 1) qgram_p[qgram] += 1;
    }
    for (uint64_t qi = r_pos, qgram = 0; qi < r_pos + r_len; qi++) {
        qgram = ((qgram << 2) | qdna(r[qi])) & MASK;
        if (qi - r_pos >= QG - 1) qgram_r[qgram] += 1;
    }
    int dist = 0;
    for (uint64_t qi = 0; qi < QSZ; qi++)  {
        dist += min(qgram_p[qi], qgram_r[qi]);
        qgram_p[qi] = qgram_r[qi] = 0;
    }

    if (dist < minqg) {
        QGRAM_NORMAL_FAILED++;
        return false;
    }

    map<int, int> hits;
    assert(q_len == r_len);
    int common = 0;
    aho->search(q.c_str() + q_pos, q_len, hits, 1);
    aho->search(r.c_str() + r_pos, r_len, hits, 2);
    for (auto &h: hits) if (h.second == 3) common++;

    double boundary = (1.0/2) * (q_len / 50.0);
    if (common < boundary) {
        CORE_FAILED++;
        return false;
    }

    return true;
}
