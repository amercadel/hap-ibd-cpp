#ifndef HAP_IBD_HPP
#define HAP_IBD_HPP
#include "match.hpp"
#include "params.hpp"



class hapIBDCpp{
	public:
		hapIBDCpp(hap_ibd_parameters params){
			this->input_vcf = params.input_vcf;
			this->plink_rate_map = params.map_file;
			this->output_file_path = params.output_file;
			if(params.use_hash_set){
				getSiteMappingAndGenotypes(this->input_vcf, this->alt_map, this->site_mapping);
			}else{
				getSiteMappingAndGenotypes(this->input_vcf, this->genotype_array, this->site_mapping);
			}			
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
			this->use_hash_set = params.use_hash_set;
			pbwtInit();
			this->p = 0;
			if(this->p){
				pbwtDestroy(this->p);
				exit(1);
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
		bool use_hash_set = false;


		std::vector<Match> matches;
		rateMapData gen_map;
		std::vector<int> site_mapping;
		std::vector<std::vector<int>> genotype_array;
		std::vector<std::unordered_set<int>> alt_map;

		std::unordered_set<std::string> output_strs;
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
					if(this->use_hash_set){
						extendBoundaryStart(m.hap1, m.hap2, m.start_site, this->alt_map);
					}
					else{
						m.start_site = extendBoundaryStart(m.hap1, m.hap2, m.start_site, this->genotype_array);
					}
					

				}
				if(m.end_site == this->windows[index].second){
					if(this->use_hash_set){
						m.end_site = extendBoundaryEnd(m.hap1, m.hap2, m.end_site, this->site_mapping, this->alt_map);
					}
					else{
						m.end_site = extendBoundaryEnd(m.hap1, m.hap2, m.end_site, this->site_mapping, this->genotype_array);
					}
					
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
				if(this->use_hash_set){
					out = processSeed(m.hap1, m.hap2, m.start_site, m.end_site, this->max_gap, this->site_mapping, this->matches, this->gen_map, this->min_seed, this->min_extend, this->min_markers, this->min_markers_extend, this->min_output, this->alt_map);
				}
				else{
					out = processSeed(m.hap1, m.hap2, m.start_site, m.end_site, this->max_gap, this->site_mapping, this->matches, this->gen_map, this->min_seed, this->min_extend, this->min_markers, this->min_markers_extend, this->min_output, this->genotype_array);
				}
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
	std::cout << "  -h use-hash-set default: false" << std::endl;
}




#endif
