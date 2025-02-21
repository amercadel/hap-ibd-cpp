#ifndef READ_RATE_MAP_HPP
#define READ_RATE_MAP_HPP
#include "utils.hpp"



struct rateMapData{
    std::vector<int> bp_vec;
    std::vector<double> cm_vec;
    std::vector<double> interpolated_cm;
    
    std::vector<double> interpolateVector(std::vector<int> &sites);
    double interpolateBasePairToGenPos(int site);
    double genPos(int site);
    int last_bp;
    double last_cm;

    void display(int idx){
        std::cout << "bp value: " << bp_vec[idx] << std::endl << "cm value: " << cm_vec[idx] << std::endl;
    };

};
rateMapData readRateMap(char* filename, std::vector<int> &sites);

// gets genetic position based on interpolated genetic distance
// OUTPUT: a double representing a site's genetic position
double getGeneticPosition(std::vector<double> &interpolated_cm, int site_idx);

// simple splitting function
// OUTPUT: a vector of strings, split based on delim
std::vector<std::string> split(std::string &line, char delim);

#endif
