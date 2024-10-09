#pragma once
#include <vector>
#include <utility>


template<typename T>
T roundToNDigits(T num, int n_digits){
    auto scale = std::pow(10.0, n_digits);
    return std::round(num * scale) / scale;
}

template<typename T>
int findVectorIndex(std::vector<T> &vec, T val){
    int left = 0;
    int right = vec.size() - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (vec[mid] == val) {
            return mid;
        } else if (vec[mid] < val) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    return -1;
}
template<typename T>
int findInsertionIndex(std::vector<T> &vec, T val){
    int left = 0;
    int right = vec.size() - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (vec[mid] == val) {
            return mid;
        } else if (vec[mid] < val) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    return left;
}
template<typename T>
int findInsertionIndex(std::vector<T> &vec, int from, int to, T val){
    int left = from;
    int right = to;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (vec[mid] == val) {
            return mid;
        } else if (vec[mid] < val) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    return left;
}


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


