#ifndef UTILS_HPP
#define UTILS_HPP
#include <vector>
#include <cmath>
#include <string>
#include <utility>
#include <sstream>
#include <fstream>
#include <iostream>
#include <chrono>
#include <thread>

extern "C"{
#include "array.h"
#include "pbwt.h"
}
    

// simple binary search for an insertion index if you need to only search a certain subset of a vector
// OUTPUT: where a value would be inserted into vector to keep it sorted
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

// simple binary search for an insertion index across and entire vector
// OUTPUT: where a value would be inserted into vector to keep it sorted
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


// simple binary search for finding a value, 
// OUTPUT: index if found, -1 if not
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

// simple rounding function
template<typename T>
T roundToNDigits(T num, int n_digits){
    auto scale = pow(10.0, n_digits);
    return round(num * scale) / scale;
}

// generates overlapping windows to guarantee that a seed is found, but allowing for parallelization
// OUTPUT: a vector of pairs representing where to split the VCF to guarantee seed will be found
std::vector<std::pair<int, int>> overlappingWindows(std::vector<double> &cm, double min_seed, int n_threads);
int minSites(std::vector<double> &cm_mapping, double min_seed);


// creates an array (using Durbin's in-house implementation in order to use for subsetting a PBWT struct)
// OUTPUT: an Array ranging from start to end (both inclusive)
Array createRangeArray(PBWT* p, int start, int end);

// simple function for outputting a bool to a string
std::string boolToString(bool b);

#endif



