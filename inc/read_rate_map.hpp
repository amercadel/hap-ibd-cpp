#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <unordered_map>

std::vector<std::string> split(std::string &line, char delim){
    std::vector<std::string> split_string;
    std::string str;
    std::stringstream ss(line);
    while (std::getline(ss, str, delim)){
        split_string.push_back(str);
    }
    return split_string;

}

template<typename T>
T findVectorIndex(std::vector<T> &vec, T val){
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
T findInsertionIndex(std::vector<T> &vec, T val){
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

struct rateMapData{
    std::vector<int> bp_vec;
    std::vector<float> cm_vec;
    std::vector<float> interpolated_cm;
    
    std::vector<float> interpolateVector(std::vector<int> &sites, std::vector<int> &bp_vec, std::vector<float> &cm_vec);
    float interpolateBasePairToGenPos(int site);
    int last_bp;
    float last_cm;

    void display(int idx){
        std::cout << "bp value: " << bp_vec[idx] << std::endl << "cm value: " << cm_vec[idx] << std::endl;
    };

};
std::vector<float> rateMapData::interpolateVector(std::vector<int> &sites, std::vector<int> &bp, std::vector<float> &cm) {
    std::vector<float> interpolated_cm;
    int n = sites.size();
    int m = bp.size();
    for (int i = 0; i < n; i++){
        int j = findInsertionIndex(bp_vec, sites[i]);
        if (j == 0){
            interpolated_cm.push_back(cm[0]);
        }else if (j == m){
            interpolated_cm.push_back(cm[m - 1]);
        }else{
            int x0 = bp[j - 1];
            int x1 = bp[j];
            float y0 = cm[j - 1];
            float y1 = cm[j];
            int x = sites[i];
            float y = ((y0 * (x1 - x)) + (y1 * (x - x0))) / (x1 - x0);
            interpolated_cm.push_back(y);
        }
    }
    return interpolated_cm;
}


rateMapData readRateMap(char* filename, std::vector<int> &sites){
    std::ifstream inputFile;
    std::string line;
    std::vector<int> bp_vec;
    std::vector<float> cm_vec;
    inputFile.open(filename);
    std::cout << "reading genetic map file\n";
    rateMapData res;
    while (getline(inputFile, line)) {
        int bp;
        float cm;
        std::string tmp;
        std::stringstream input_str(line); 
        std::vector<std::string> data = split(line, ' ');
        bp = std::stoi(data[3]);
        cm = std::stof(data[2]);
        bp_vec.push_back(bp);
        cm_vec.push_back(cm);

    }
    std::cout << "rate map loaded into memory\n";
    res.bp_vec = bp_vec;
    res.cm_vec = cm_vec;
    res.last_bp = bp_vec[bp_vec.size() - 1];
    res.last_cm = cm_vec[cm_vec.size() - 1];
    std::cout << "interpolating rate map\n";
    res.interpolated_cm = res.interpolateVector(sites, res.bp_vec, res.cm_vec);
    std::cout << "rate map interpolated\n";
    return res;
}
