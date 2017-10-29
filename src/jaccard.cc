/// 786 


/// start Thu Oct 26 10:56:02 PDT 2017
/// 

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "common.h"
#include "search.h"
using namespace std;

extern int QGRAM_NORMAL_FAILED;
extern int QGRAM_SPACED_FAILED;
extern int CORE_FAILED;

int main(int argc, char **argv)
{
    if (argc < 3) exit(1);

    eprn("🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁", 0);
    eprn("🐚       S   E   D   E   F       🐚", 0);
    eprn("🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁", 0);
    eprnn("   🖥  {}; arguments: ", GITVER);
    for (int i = 0; i < argc; i++) eprnn(" {}", argv[i]);
    eprnn("\n");
    
    string ref_path = argv[1];
    string query_path = argv[2];
    bool is_complement = (argc > 3 && tolower(argv[3][0]) == 'y');

    string line, reference, query;

    ifstream fin;

    fin.open(ref_path.c_str());
    string ref_chr;
    while (getline(fin, line)) {
        if (line[0] != '>') reference += line;
        else ref_chr = line.substr(1);
    }
    fin.close();
    fin.clear();

    transform(reference.begin(), reference.end(), reference.begin(), ::toupper);
    if (is_complement) {
        eprnn("Reversing reference...\n");
        reverse(reference.begin(), reference.end());
        for (int i = 0; i < reference.size(); i++) {
            char &c = reference[i];
            c = (c == 'A' ? 'T' : (c == 'C' ? 'G' : (c == 'G' ? 'C' : (c == 'T' ? 'A' : c))));
        }
    }

    fin.open(query_path.c_str());
    string query_chr;
    while (getline(fin, line)) {
        if (line[0] != '>') query += line;
        else query_chr = line.substr(1);
    }
    fin.close();
    transform(query.begin(), query.end(), query.begin(), ::toupper);

    Hash ref_hash = Hash(reference);

    Hash *query_hash = &ref_hash; 
    Hash tmp;
    if (query_path != ref_path || query.size() != reference.size() || is_complement) {
        tmp = Hash(query);
        query_hash = &tmp;
    }

    bool allow_overlaps = (ref_path != query_path) || is_complement;
    eprn("Allowing overlaps: {}", allow_overlaps);
    eprn("Reverse complement: {}", is_complement);

    int total = 0;
    for (int i = 0; i < query.size(); i += 250) {
    // chr22    16467109    16469072    chr22:16883317  0   _   chr22   16883317    16885246
    // 34421249    34419320
    // for (int i = 16467100, j = 0; i < query.size(); i += 250, j++) {
        while (query[i] == 'N') i++;
        while (i % 250 != 0) i++;
        if (i % 5000 == 0) {
            double perc = 100.0 * i / double(query.size());
            eprnn("\r   {} {:.1f}% ({})", string(int(perc / 2) + 1, '-'), perc, i);
        }
        // prn("{}", i);

        vector<Hit> mapping = search(i, ref_hash, *query_hash, allow_overlaps);
        for (int i = 0; i < mapping.size(); i++) {
            Hit &pp = mapping[i];
            // BEDPE
            if (is_complement) {
                pp.i = ref_hash.seq.size() - pp.i + 1;
                pp.j = ref_hash.seq.size() - pp.j + 1;
                swap(pp.i, pp.j);
            }
            prn("{}\t{}\t{}\t{}\t{}\t{}\t\t{:.0f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                query_chr, pp.p, pp.q, 
                ref_chr, pp.i, pp.j,
                pp.id, 
                "+", is_complement ? "-" : "+",
                // Optional fields
                int(pp.break_criteria),
                pp.q - pp.p ,
                pp.jaccard.first, pp.jaccard.second,
                pp.edist.first, pp.edist.second
            );
            
            total += 1;
        }
    }
    eprnn("\n");

    eprn("Total:               {:10n}", total);
    eprn("Fails: cores         {:10n}\n"
         "       q-gram normal {:10n}\n"
         "       q-gram spaced {:10n}", CORE_FAILED, QGRAM_NORMAL_FAILED, QGRAM_SPACED_FAILED);

    return 0;
}