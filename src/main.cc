/// 786 

/******************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <locale>

#include "jaccard.h"
#include "aho.h"
#include "align.h"
#include "wgac.h"

using namespace std;

/******************************************************************************/

extern shared_ptr<AHOAutomata> aho;

/******************************************************************************/

void parse(int argc, char **argv);

/******************************************************************************/

int main(int argc, char **argv)
{
	ios_base::sync_with_stdio(0);
	setlocale(LC_NUMERIC, "en_US.UTF-8");
	if (argc < 3) exit(1);

	eprn("🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁");
	eprn("🐚       S   E   D   E   F       🐚");
	eprn("🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁 🍁");
	eprnn("   🖥  {}; arguments: ", GITVER);
	for (int i = 0; i < argc; i++) eprnn(" {}", argv[i]);
	eprn("");
	
	string command = argv[1];

	aho = make_shared<AHOAutomata>();

	if (command == "search") {
		string ref_path = argv[2];
		string query_path = argv[3];
		bool is_complement = (argc > 4 && tolower(argv[4][0]) == 'y');
		jaccard_search(ref_path, query_path, is_complement);
	} else if (command == "align") {
		string ref_path = argv[2];
		string bed_path = argv[3];
		int resume_after = argc <= 4 ? -1 : atoi(argv[4]);
		align_main(ref_path, bed_path, resume_after);
	} else if (command == "parse") {
		parse(argc-1, argv+1);
	} else if (command == "wgac") {
		check_wgac(argv[2]);
	} else {
		eprn("Whoops, invalid command!");
	}
	
	return 0;
}