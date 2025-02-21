#include "hap_ibd.hpp"

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
		} else if (std::string(argv[i]) == "-h") {
			params.use_hash_set = true;
		}
	}

	displayParameters(params);
	hapIBDCpp obj(params);
	
	auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
	
    std::cout << "time elapsed: " << std::chrono::duration<double>(diff).count() << " seconds" << std::endl;
	



	
	
	return 0;

}


