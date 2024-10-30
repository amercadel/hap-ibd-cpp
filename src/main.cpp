#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include <set>
#include <filesystem>
#include "match.hpp"
#include "vcf.hpp"
#include "utils.hpp"
#include "omp.h"
extern "C"
{
#include "pbwt.h"
}
// need to fix issue with min markers extend vs min markers seed. first I'd like to clean everything up

class hapIBDCpp{
	public:
		hapIBDCpp(char* input_vcf, char* plink_rate_map, const char* output_file_path = "output.txt", double min_seed = 2.0f, int max_gap = 1000, 
					double min_extend = 1.0f, double min_output = 2.0f, int min_markers = 100, int min_mac = 2, int n_threads = 1){
			this->input_vcf = input_vcf;
			this->plink_rate_map = plink_rate_map;
			this->output_file_path = output_file_path;
			getSiteMappingAndGenotypes(input_vcf, this->genotype_array, this->site_mapping);
			this->gen_map = readRateMap(plink_rate_map, site_mapping);
			this->min_seed = min_seed;
			this->max_gap = max_gap;
			this->min_extend = std::min(1.0, this->min_seed);
			this->min_output = min_output;
			this->min_markers = min_markers;
			this->min_mac = min_mac;
			this->min_markers_extend = floor((this->min_extend/this->min_seed) * this->min_markers);
			this->n_threads = n_threads;


			this->windows = overlappingWindows(this->gen_map.interpolated_cm, this->min_seed, this->min_markers, this->n_threads);
			this->intermediate_files = splitVCFByPos(this->input_vcf, windows);

			std::vector<std::thread> threads;
			for(int i = 0; i < this->n_threads; ++i){
				std::thread t([this, i]() { run(intermediate_files[i], i); });
				threads.push_back(std::move(t)); // Use std::move to move the thread object
			}
			for(auto& t: threads){
				t.join();
			}
			outputSegments();
	
	
			
	}



	private:
		char* input_vcf;
		char* plink_rate_map;
		const char* output_file_path;
		double min_seed;
		int max_gap;
		double min_extend;
		double min_output;
		int min_markers;
		int min_mac;
		int min_markers_extend;
		int n_threads;
		std::vector<char*> intermediate_files;
		std::vector<Match> matches;
		rateMapData gen_map;
		std::vector<int> site_mapping;
		std::vector<std::vector<int>> genotype_array;
		std::set<std::string> output_strs;
		std::vector<std::pair<int, int>> windows;
		
		int* runPBWT(char* input_vcf, int index){
			pbwtInit();
			PBWT* p = 0;
			if(p){
				pbwtDestroy(p);
			}
			p = pbwtReadVcfGT(input_vcf);
			int* raw_matches = pbwtLongMatches(p, 100, index);
			return raw_matches;

			

		}

		std::vector<Match> getMatches(int* matches_array, int index){
			int i = 0;
			std::vector<Match> seeds;
			while(matches_array[i] != -1){
				Match m(matches_array[i], matches_array[i+1], matches_array[i+2] + this->windows[index].first, matches_array[i+3] + this->windows[index].first);
				double f1 = getGeneticPosition(gen_map.interpolated_cm, m.start_site);
				double f2 = getGeneticPosition(gen_map.interpolated_cm, m.end_site);
				double len = f2 - f1;
				if ((len >= this->min_seed) && ((m.n_sites) >= this->min_markers)){
					m.len_cm = len;
					seeds.push_back(m);
					}
				i = i + 4;
			}
			return seeds;

		}
		void processSeeds(std::vector<Match> seeds){
			for(size_t c = 0; c < seeds.size(); c++){				
				
				Match m = seeds[c];
				std::string out;
				out = processSeed(m.hap1, m.hap2, m.start_site, m.end_site, this->max_gap, this->site_mapping, this->matches, this->gen_map, this->min_seed, this->min_extend, this->min_markers, this->min_markers_extend, this->min_output, this->genotype_array);
				if(!out.empty()){
					this->output_strs.insert(out);
				}
			}
			
		}

		void run(char* split_vcf, int index){
			int* raw_matches = runPBWT(this->intermediate_files[index], index);
			std::vector<Match> filtered_matches = getMatches(raw_matches, index);
			processSeeds(filtered_matches);

		}

		void outputSegments(){
			std::ofstream output_file(this->output_file_path);
			if (output_file.is_open()) {
				for (const auto& str : this->output_strs) {
					output_file << str;
				}
				output_file.close();
			} else {
				std::cerr << "Unable to open file: " << this->output_file_path << std::endl;
			}
			for(int i = 0; i < this->n_threads; i++){
				bool ret = std::filesystem::remove(this->intermediate_files[i]);
				if(!ret){
					std::cerr << "Error removing intermediate files\n";
				}
			}
		}
};



int main(int argc, char **argv){
	auto start = std::chrono::steady_clock::now();
	char *input_vcf = argv[1];
	char *plink_rate_map = argv[2];
	
	hapIBDCpp obj(input_vcf, plink_rate_map);
	
	
	auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
	
    std::cout << std::chrono::duration<double>(diff).count() << " seconds" << std::endl;

	
	
	return 0;

}


