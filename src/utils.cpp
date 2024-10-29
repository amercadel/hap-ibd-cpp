#include <utility>
#include <string>
#include <sstream>
#include "utils.hpp"

std::vector<std::pair<int, int>> overlappingWindows(std::vector<double> cm, double min_seed, int min_markers, int n_threads){
    std::vector<std::pair<int, int>> window_vec;
    double L = cm.back() - cm.front();
    double window_size = ((L - min_seed) / n_threads + min_seed);
    double dist = (L - min_seed) / n_threads;
    double start_gen_pos = 0.0f;
    int i = 0;
    while(i < n_threads){
        double end_gen_pos = start_gen_pos + window_size;
        int start = findInsertionIndex(cm, start_gen_pos);
        int end = findInsertionIndex(cm, end_gen_pos);
        start_gen_pos += dist;
        window_vec.push_back(std::pair<int, int>(start, end));
        i++;
    }
    return window_vec;
}


std::vector<std::string> getIntermediateMatchFileNames(int n_threads){
    std::vector<std::string> file_names;
    std::string name;
    for(int i = 0; i < n_threads; i++){
        std::stringstream ss;
        ss << "intermediate_matches_" << i << ".txt";
        name = ss.str();
        file_names.push_back(name);
    }
    return file_names;
}

class HapIBDParameters{
    // 10/29: will implement later
    public:
        const char* input_vcf;
        const char* plink_rate_map;
        const char* output_file_path;
        double min_seed;
        int max_gap;
        double min_extend;
        double min_output;
        int min_markers;
        int n_threads;
    private:
        int min_markers_extend;


};

