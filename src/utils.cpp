#include <utility>
#include <string>
#include <sstream>
#include "utils.hpp"

std::vector<std::pair<int, int>> overlappingWindows(std::vector<double> &cm, double min_seed, int n_threads){
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
int minSites(std::vector<double> &cm_mapping, double min_seed) {
    int min = cm_mapping.size();
    int j = 100; // heuristic: it is unlikely that min seed of 2 would cover just 2 sites or even 10. Maybe start at 50? or 100?
    for (int i = 0; i < cm_mapping.size(); i++) {
        while (j < cm_mapping.size() && cm_mapping[j] - cm_mapping[i] < min_seed) {
            j = j + 2;
        }
        if (cm_mapping.back() - cm_mapping[i] > min_seed) {
            int n = j - i;
            if (n < min) {
                min = n;
            }
        }
    }
    return (min == cm_mapping.size()) ? -1 : min;
}


