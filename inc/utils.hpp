#pragma once
#include <vector>

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
T roundToNDigits(T num, int n_digits){
    auto scale = std::pow(10.0, n_digits);
    return std::round(num * scale) / scale;
}


std::vector<std::pair<int, int>> overlappingWindows(std::vector<float> cm, float min_seed, int min_markers, int n_threads);


