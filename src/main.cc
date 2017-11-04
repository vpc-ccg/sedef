/// 786 

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "jaccard.h"
#include "align.h"

using namespace std;

void parse(int argc, char **argv);

int main(int argc, char **argv)
{
	std::ios_base::sync_with_stdio(0);
    //if (argc < 3) exit(1);

    eprn("🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁", 0);
    eprn("🐚       S   E   D   E   F       🐚", 0);
    eprn("🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁", 0);
    eprnn("   🖥  {}; arguments: ", GITVER);
    for (int i = 0; i < argc; i++) eprnn(" {}", argv[i]);
    eprnn("\n");
    
	string command = argv[1];

	if (command == "search") {
    	string ref_path = argv[2];
    	string query_path = argv[3];
    	bool is_complement = (argc > 4 && tolower(argv[4][0]) == 'y');
    	jaccard_search(ref_path, query_path, is_complement);
    } else if (command == "align") {
    	string ref_path = argv[2];
		string bed_path = argv[3];
    	int resume_after = argc <= 4 ? -1 : atoi(argv[4]);
    	align(ref_path, bed_path, resume_after);
    } else if (command == "parse") {
        parse(argc-1, argv+1);
    } else {
    	eprnn("Whoops, invalid command!\n");
    }
    
    return 0;
}