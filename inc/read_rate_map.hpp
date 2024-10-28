#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cassert>
#include "utils.hpp"



struct rateMapData{
    std::vector<int> bp_vec;
    std::vector<float> cm_vec;
    std::vector<float> interpolated_cm;
    
    std::vector<float> interpolateVector(std::vector<int> &sites);
    float interpolateBasePairToGenPos(int site);
    float genPos(int site);
    int last_bp;
    float last_cm;

    void display(int idx){
        std::cout << "bp value: " << bp_vec[idx] << std::endl << "cm value: " << cm_vec[idx] << std::endl;
    };

};
rateMapData readRateMap(char* filename, std::vector<int> &sites);
float getGeneticPosition(std::vector<float> &interpolated_cm, int site_idx);
std::vector<std::string> split(std::string &line, char delim);
