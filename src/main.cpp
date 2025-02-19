#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <set>
#include <filesystem>
#include <unordered_set>

#include "match.hpp"
#include "vcf.hpp"
#include "utils.hpp"
#include "params.hpp"
extern "C"
{
#include "array.h"
}


class hapIBDCpp{
	public:
		hapIBDCpp(hap_ibd_parameters params){
			this->input_vcf = params.input_vcf;
			this->plink_rate_map = params.map_file;
			this->output_file_path = params.output_file;
			getSiteMappingAndGenotypes(this->input_vcf, this->genotype_array, this->site_mapping, this->n_threads);
			// getSiteMappingAndGenotypes(this->input_vcf, this->alt_map, this->site_mapping, this->n_threads);
			this->gen_map = readRateMap(this->plink_rate_map, site_mapping);
			this->min_seed = params.min_seed;
			this->max_gap = params.max_gap;
			this->min_extend = std::min(1.0, params.min_seed);
			this->min_output = params.min_output;
			this->min_markers = params.min_markers;
			this->min_mac = params.min_mac;
			this->min_markers_extend = floor((this->min_extend/this->min_seed) * this->min_markers);
			this->n_threads = params.n_threads;
			this->cm_threshold_sites = minSites(this->gen_map.interpolated_cm, this->min_seed);
			pbwtInit();
			this->p = 0;
			if(this->p){
				pbwtDestroy(this->p);
				die("Critical error");
			}
			this->p = pbwtReadVcfGT(this->input_vcf);


			this->windows = overlappingWindows(this->gen_map.interpolated_cm, this->min_seed, this->n_threads);
			


			if(this->n_threads > 1){
				

				std::vector<PBWT*> pbwt_vec;

				for(int i = 0; i < this->n_threads; i++){
					std::pair<int, int> pa = windows[i];
					Array a = createRangeArray(p, pa.first, pa.second);
					PBWT* pt = 0;
					pt = pbwtSelectSites(p, a, true);
					pbwt_vec.push_back(pt);
				}
				std::vector<std::thread> threads;
				for(int j = 0; j < this->n_threads; ++j){
					std::thread t(&hapIBDCpp::run, this, pbwt_vec[j], j);
					threads.push_back(std::move(t)); // Use std::move to move the thread object
				}
				for(auto& t: threads){
					t.join();
				}
			}
			else{
				run(p, 0);
			}
			

			outputSegments();
			
	}
	private:

		PBWT* p;

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
		int cm_threshold_sites;
		std::vector<char*> intermediate_files;
		std::vector<Match> matches;
		rateMapData gen_map;
		std::vector<int> site_mapping;
		std::vector<std::vector<int>> genotype_array;
		std::unordered_map<int, std::bitset<MAX_N_SAMPLES>> alt_map;
		std::set<std::string> output_strs;
		std::vector<std::pair<int, int>> windows;
		
		int* runPBWT(PBWT* split_p, int index){
			int* raw_matches = pbwtLongMatches(split_p, this->cm_threshold_sites, index);
			return raw_matches;

		}

		std::vector<Match> getMatches(int* matches_array, int index){
			int i = 0;
			std::vector<Match> seeds;
			while(matches_array[i] != -1){
				Match m(matches_array[i], matches_array[i+1], matches_array[i+2] + this->windows[index].first, matches_array[i+3] + this->windows[index].first);
				if((m.start_site == this->windows[index].first) && (m.start_site != 0)){ // no need to check for extension if it starts at 0
					m.start_site = extendBoundaryStart(m.hap1, m.hap2, m.start_site, this->genotype_array);
					// m.start_site = extendBoundaryStart(m.hap1, m.hap2, m.start_site, this->alt_map);

				}
				if(m.end_site == this->windows[index].second){
					m.end_site = extendBoundaryEnd(m.hap1, m.hap2, m.end_site, this->site_mapping, this->genotype_array);
					// m.end_site = extendBoundaryEnd(m.hap1, m.hap2, m.end_site, this->site_mapping, this->alt_map);
				}
				m.n_sites = m.end_site - m.start_site;
				
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
		void processSeeds(std::vector<Match> &seeds){
			for(size_t c = 0; c < seeds.size(); c++){				
				
				Match m = seeds[c];
				std::string out;
				out = processSeed(m.hap1, m.hap2, m.start_site, m.end_site, this->max_gap, this->site_mapping, this->matches, this->gen_map, this->min_seed, this->min_extend, this->min_markers, this->min_markers_extend, this->min_output, this->genotype_array);
				// out = processSeed(m.hap1, m.hap2, m.start_site, m.end_site, this->max_gap, this->site_mapping, this->matches, this->gen_map, this->min_seed, this->min_extend, this->min_markers, this->min_markers_extend, this->min_output, this->alt_map);
				if(!out.empty()){
					this->output_strs.insert(out);
				}
				
			}
			
		}

		void run(PBWT* split_p, int index){
			int* raw_matches = runPBWT(split_p, index);
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
		}
};


void printUsage(){
	std::cout << "hap-ibd-cpp\n";
	std::cout << "Usage: ./hap-ibd -i {input-vcf} -m {plink-rate-map-file} -o {output-file-path}" << std::endl;
	std::cout << "Optional Parameters:" << std::endl;
	std::cout << "  -s: min-seed    default: 2.0" << std::endl;
	std::cout << "  -g: max-gap     default: 1000" << std::endl;
	std::cout << "  -e min-extend   default: 1.0" << std::endl;
	std::cout << "  -m min-output   default: 2.0" << std::endl;
	std::cout << "  -k min-markers  default: 100" << std::endl;
	std::cout << "  -a min-mac      default: 2" << std::endl;
	std::cout << "  -t n-threads    default: 1" << std::endl;
}


int main(int argc, char **argv){
	
	auto start = std::chrono::steady_clock::now();
	
	if (argc == 1){
		printUsage();
		return 0;
	}
	hap_ibd_parameters params;
	for(int i = 0; i < argc; i++){
		if (std::string(argv[i]) == "-i") {
			params.input_vcf = argv[++i];
		} else if (std::string(argv[i]) == "-m") {
			params.map_file = argv[++i];
		} else if (std::string(argv[i]) == "-o") {
			params.output_file = argv[++i];
		} else if (std::string(argv[i]) == "-s") {
			params.min_seed = std::stod(argv[++i]);
		} else if (std::string(argv[i]) == "-g") {
			params.max_gap = std::stoi(argv[++i]);
		} else if (std::string(argv[i]) == "-e") {
			params.min_extend = std::stod(argv[++i]);
		} else if (std::string(argv[i]) == "-m") {
			params.min_output = std::stod(argv[++i]);
		} else if (std::string(argv[i]) == "-k") {
			params.min_markers = std::stoi(argv[++i]);
		} else if (std::string(argv[i]) == "-a") {
			params.min_mac = std::stoi(argv[++i]);
		} else if (std::string(argv[i]) == "-t") {
			params.n_threads = std::stoi(argv[++i]);
		}
	}

	std::cout << "input-vcf: " << params.input_vcf << std::endl;
	std::cout << "plink-rate-map-file: " << params.map_file << std::endl;
	std::cout << "output-file-path: " << params.output_file << std::endl;
	std::cout << "min-seed: " << params.min_seed << std::endl;
	std::cout << "max-gap: " << params.max_gap << std::endl;
	std::cout << "min-extend: " << params.min_extend << std::endl;
	std::cout << "min-output: " << params.min_output << std::endl;
	std::cout << "min-markers: " << params.min_markers << std::endl;
	std::cout << "min-mac: " << params.min_mac << std::endl;
	std::cout << "n-threads: " << params.n_threads << std::endl;
	
	hapIBDCpp obj(params);
	
	
	auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
	
    std::cout << "time elapsed: " << std::chrono::duration<double>(diff).count() << " seconds" << std::endl;
	



	
	
	return 0;

}


