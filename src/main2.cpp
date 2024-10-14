#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <unistd.h>
#include "match.hpp"
#include "vcf.hpp"
extern "C"
{
#include "pbwt.h"
}


class hapIBDCpp{
	public:
		hapIBDCpp(char* input_vcf, char* plink_rate_map, char* output_file_path = "output.txt", float min_seed = 2.0, int max_gap = 100, 
					float min_extend = 1.0, float min_output = 2.0, int min_markers = 100, int min_mac = 2){
			this->input_vcf = input_vcf;
			this->plink_rate_map = plink_rate_map;
			this->output_file_path = output_file_path;
			this->site_mapping = getSiteMapping(input_vcf);
			this->gen_map = readRateMap(plink_rate_map, site_mapping);
			this->min_seed = min_seed;
			this->max_gap = max_gap;
			this->min_extend = std::min(1.0f, this->min_seed);
			this->min_output = min_output;
			this->min_markers = min_markers;
			this->min_mac = min_mac;
			
			std::vector<std::pair<int, int>> windows = overlappingWindows(this->gen_map.interpolated_cm, this->min_seed, this->min_markers, 4);
			std::vector<char*> intermediate_files = splitVCFByPos(this->input_vcf, windows);

			
			for(size_t i = 0; i < intermediate_files.size(); i++){
				runPBWT(intermediate_files[i], i);
			}
			// std::cout << "running pbwt\n";
			// this->runPBWT(this->input_vcf);
			std::cout << "fetching matches" << std::endl;
			this->getMatches();
			std::cout << "processing seed\n";
			this->processSeeds();
			
			
			
			

			
		}



	private:
		char* input_vcf;
		char* plink_rate_map;
		char* output_file_path;
		char* match_file = "intermediate_matches.txt";
		float min_seed;
		int max_gap;
		float min_extend;
		float min_output;
		int min_markers;
		int min_mac;
		std::vector<Match> filtered_matches;
		rateMapData gen_map;
		std::vector<int> site_mapping;
		void runPBWT(char* input_vcf, int index){
			pbwtInit();
			PBWT* p = 0;
			if(p){
				pbwtDestroy(p);
			}
			p = pbwtReadVcfGT(input_vcf);
			pbwtLongMatches(p, 0, index);

			

		}
		void getMatches(){
			std::vector<Match> *matches = new std::vector<Match>;
			std::ifstream mf;
			mf.open(this->match_file);
			std::string line;
			
			while(std::getline(mf, line)){
				Match m(line);
				
				
				// std::cout << m.start_site << " " << m.end_site << std::endl;
				float f1 = getGeneticPosition(gen_map.interpolated_cm, m.start_site);
				float f2 = getGeneticPosition(gen_map.interpolated_cm, m.end_site - 1);
				

				float len = f2 - f1;
				if (len >= this->min_extend){
					m.len_cm = len;
					
					matches->push_back(m);
					}
			}
			mf.close();
			std::cout << "got matches\n";
			size_t i = 0;
			while(i < matches->size() - 1){
				if((*matches)[i] == (*matches)[i+1]){
					this->filtered_matches.push_back((*matches)[i]);
					i = i + 2;
				}
				else{
					this->filtered_matches.push_back((*matches)[i]);
					i++;
				}
			}
			std::sort(filtered_matches.begin(), filtered_matches.end(), compareHaps);
			delete matches;
		}
		void processSeeds(){
			std::ofstream output_file(this->output_file_path);
			for(size_t c = 0; c < filtered_matches.size(); c++){
				Match m = filtered_matches[c];
				std::string out;
				out = processSeed(m.hap1, m.hap2, m.start_site, m.end_site, this->max_gap, this->site_mapping, this->filtered_matches, this->gen_map, this->min_seed, this->min_extend, this->min_markers, this->min_output);
				output_file << out;
			}
			output_file.close();


		}
};



int main(int argc, char **argv){
	auto start = std::chrono::steady_clock::now();
	char *input_vcf = argv[1];
	char *plink_rate_map = argv[2];
	
	hapIBDCpp run(input_vcf, plink_rate_map);
	auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << std::chrono::duration<double>(diff).count() << " seconds" << std::endl;

	
	
	return 0;

}


