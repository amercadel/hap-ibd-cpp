#ifndef PARAMS_HPP
#define PARAMS_HPP
#include "utils.hpp"

typedef struct hap_ibd_parameters {
    char* input_vcf;
    char* map_file;
    char* output_file;
    double min_seed = 2.0;
    int max_gap = 1000;
    double min_extend = 1.0;
    double min_output = 2.0;
    int min_mac = 2;
    int min_markers = 100;
    int n_threads = 1;
    bool use_hash_set = false;
} param_t;


void displayParameters(param_t params){
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
	std::cout << "using-hash-set: " << boolToString(params.use_hash_set) << std::endl;
}

#endif