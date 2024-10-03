#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <unordered_map>
#include "utils.hpp"

float getGeneticPosition(std::vector<float> &interpolated_cm, int site_idx){
	float genetic_position = interpolated_cm[site_idx];
	return genetic_position;
}

std::vector<std::string> split(std::string &line, char delim){
    std::vector<std::string> split_string;
    std::string str;
    std::stringstream ss(line);
    while (std::getline(ss, str, delim)){
        split_string.push_back(str);
    }
    return split_string;

}



struct rateMapData{
    std::vector<int> bp_vec;
    std::vector<float> cm_vec;
    std::vector<float> interpolated_cm;
    
    std::vector<float> interpolateVector(std::vector<int> &sites);
    float interpolateBasePairToGenPos(int site);
    int last_bp;
    float last_cm;

    void display(int idx){
        std::cout << "bp value: " << bp_vec[idx] << std::endl << "cm value: " << cm_vec[idx] << std::endl;
    };

};
float rateMapData::interpolateBasePairToGenPos(int site){
    int idx = findVectorIndex(bp_vec, site);
    if (idx > 0){
        return cm_vec[idx];
    }
    else{
        idx = findInsertionIndex(bp_vec, site);
        if (idx == 0){
            return cm_vec[0];
        }
        else if(idx == bp_vec.size()){
            return cm_vec.back();
        }
        else{
            int x0 = bp_vec[idx - 1];
            int x1 = bp_vec[idx];
            float y0 = cm_vec[idx - 1];
            float y1 = cm_vec[idx];
            float y = ((y0 * (x1 - site)) + (y1 * (site - x0))) / (x1 - x0);
            return y;

        }
    }


}

std::vector<float> rateMapData::interpolateVector(std::vector<int> &sites) {
    std::vector<float> interpolated_cm;
    for(size_t c = 0; c < sites.size(); c++){
        interpolated_cm.push_back(interpolateBasePairToGenPos(sites[c]));
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
    res.interpolated_cm = res.interpolateVector(sites);
    std::cout << "rate map interpolated\n";
    return res;
}
