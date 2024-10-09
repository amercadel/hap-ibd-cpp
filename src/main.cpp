#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <unistd.h>
#include "match.hpp"
#include "read_rate_map.hpp"
#include "vcf.hpp"
extern "C"
{
#include "pbwt.h"
}





int main(int argc, char **argv){
	auto start = std::chrono::steady_clock::now();
	char *input_vcf = argv[1];
	char *plink_rate_map = argv[2];
	std::vector<int> site_mapping = getSiteMapping(input_vcf);
	

	rateMapData gen_map = readRateMap(plink_rate_map, site_mapping);
	



	pbwtInit();
	PBWT* p = 0;
	if(p){
		pbwtDestroy(p);
	}
	std::cout << "running pbwt\n";
	p = pbwtReadVcfGT(input_vcf);

	
	pbwtLongMatches(p, 0);


	float min_extend = 1.0f;
    std::ifstream match_file;
	match_file.open("intermediate_matches_native.txt");
	std::string line;
	std::vector<Match> matches;
	if (!match_file.is_open()) {
        std::cerr << "Could not open the file!" << std::endl;
        return 1;
    }
	std::cout << "fetching matches\n";
	while(std::getline(match_file, line)){
		Match m(line);
		float f1 = getGeneticPosition(gen_map.interpolated_cm, m.start_site);
		float f2 = getGeneticPosition(gen_map.interpolated_cm, m.end_site);
		float len = f2 - f1;
		if (len > min_extend){
			m.len_cm = len;
			matches.push_back(m);
		}
	}
	
	match_file.close();
	std::vector<Match> filtered_matches;
	size_t i = 0;
	while(i < matches.size() - 1){
		if(matches[i] == matches[i+1]){
			filtered_matches.push_back(matches[i]);
			i = i + 2;
		}
		else{
			filtered_matches.push_back(matches[i]);
			i++;
		}
	}
	
	std::ofstream output_file("output_matches.txt");
	std::sort(filtered_matches.begin(), filtered_matches.end(), compareStartSite);
	for(size_t c = 0; c < filtered_matches.size(); c++){
		Match m = filtered_matches[c];
		std::string out;
		out = processSeed(m.hap1, m.hap2, m.start_site, m.end_site, 1000, site_mapping, filtered_matches, gen_map, 2.0f, 1.0f, 100, 2.0f);
		output_file << out;
	}
	output_file.close();
	auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << std::chrono::duration<double>(diff).count() << " seconds" << std::endl;
	
	
	return 0;

}


