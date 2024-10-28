#include <utility>
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


