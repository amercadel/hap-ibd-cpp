#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <thread>
#include <cstdio>
#include <unistd.h>
#include <unordered_set>
#include "match.hpp"
#include "vcf.hpp"
#include "utils.hpp"
// #include "omp.h"
extern "C"
{
#include "pbwt.h"
}
// need to fix issue with min markers extend vs min markers seed. first I'd like to clean everything up

class hapIBDCpp{
	public:
		hapIBDCpp(char* input_vcf, char* plink_rate_map, const char* output_file_path = "output.txt", float min_seed = 2.0f, int max_gap = 1000, 
					float min_extend = 1.0f, float min_output = 2.0f, int min_markers = 100, int min_mac = 2, int n_threads = 1){
			this->input_vcf = input_vcf;
			this->plink_rate_map = plink_rate_map;
			this->output_file_path = output_file_path;
			getSiteMappingAndGenotypes(input_vcf, this->genotype_array, this->site_mapping);
			this->gen_map = readRateMap(plink_rate_map, site_mapping);
			this->min_seed = min_seed;
			this->max_gap = max_gap;
			this->min_extend = std::min(1.0f, this->min_seed);
			this->min_output = min_output;
			this->min_markers = min_markers;
			this->min_mac = min_mac;
			this->min_markers_extend = floor((this->min_extend/this->min_seed) * this->min_markers);
			this->n_threads = n_threads;

			
			// std::vector<std::pair<int, int>> windows = overlappingWindows(this->gen_map.interpolated_cm, this->min_seed, this->min_markers, this->n_threads);
			// std::vector<char*> intermediate_files = splitVCFByPos(this->input_vcf, windows);

			
			// #pragma omp parallel for
			// for(int i = 0; i < intermediate_files.size(); i++){
			// 	runPBWT(intermediate_files[i], i);
			// }
			std::cout << "running pbwt\n";
			runPBWT(this->input_vcf, 0);
			
			std::cout << "fetching matches" << std::endl;
			this->getMatches("intermediate_matches_0.txt");
			std::cout << "processing seeds\n";
			processSeeds();
 			
			
			
	

			
			

			
		}



	private:
		char* input_vcf;
		char* plink_rate_map;
		const char* output_file_path;
		float min_seed;
		int max_gap;
		float min_extend;
		float min_output;
		int min_markers;
		int min_mac;
		int min_markers_extend;
		int n_threads;
		std::vector<Match> matches;
		rateMapData gen_map;
		std::vector<int> site_mapping;
		std::vector<std::vector<int>> genotype_array;
		void runPBWT(char* input_vcf, int index){
			pbwtInit();
			PBWT* p = 0;
			if(p){
				pbwtDestroy(p);
			}
			p = pbwtReadVcfGT(input_vcf);
			pbwtLongMatches(p, 350, index);

			

		}
		void getMatches(char* match_file){
			std::ifstream mf;
			mf.open(match_file);
			std::string line;
			while(std::getline(mf, line)){
				Match m(line);
				float f1 = getGeneticPosition(gen_map.interpolated_cm, m.start_site);
				float f2 = getGeneticPosition(gen_map.interpolated_cm, m.end_site);
				float len = f2 - f1;
				if ((len >= this->min_seed) && ((m.n_sites) >= this->min_markers)){
					m.len_cm = len;
					this->matches.push_back(m);
					}
				
			}

			std::cout << "got matches\n";

		}
		void processSeeds(){
			std::ofstream output_file(this->output_file_path);
			for(size_t c = 0; c < this->matches.size(); c++){				
				
				Match m = this->matches[c];
				std::string out;
				out = processSeed(m.hap1, m.hap2, m.start_site, m.end_site, this->max_gap, this->site_mapping, this->matches, this->gen_map, this->min_seed, this->min_extend, this->min_markers, this->min_markers_extend, this->min_output, this->genotype_array);
				if(!out.empty()){
					output_file << out;
					}
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


