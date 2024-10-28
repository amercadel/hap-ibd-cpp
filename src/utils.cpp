#include <utility>
#include <string>
#include <sstream>
#include "utils.hpp"

std::vector<std::pair<int, int>> overlappingWindows(std::vector<float> cm, float min_seed, int min_markers, int n_threads){
    std::vector<std::pair<int, int>> window_vec;
    float L = cm.back() - cm.front();
    float window_size = ((L - min_seed) / n_threads + min_seed);
    float dist = (L - min_seed) / n_threads;
    float start_gen_pos = 0.0f;
    int i = 0;
    while(i < n_threads){
        float end_gen_pos = start_gen_pos + window_size;
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

