#pragma once
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <string>
#include "read_rate_map.hpp"
#include "utils.hpp"
 




class Match{
    public:
        int hap1;
        int hap2;
        int start_site;
        int end_site;
        int n_sites;
        float len_cm;
        bool visited = false;
        Match(std::string &input_str){
            std::vector<std::string> split_str = split(input_str, '\t');
            this->hap1 = std::min(stoi(split_str[1]), stoi(split_str[2]));
            this->hap2 = std::max(stoi(split_str[1]), stoi(split_str[2]));
            this->start_site = stoi(split_str[3]);
            this->end_site = stoi(split_str[4]) - 1;
            this->n_sites = stoi(split_str[5]);
        }
        Match(int hap1, int hap2, int start_site, int end_site){
            this->hap1 = hap1;
            this->hap2 = hap2;
            this->start_site = start_site;
            this->end_site = end_site - 1;
            this->n_sites = end_site - start_site;
        }
        void display(){
            std::cout << this->hap1 << "\t" << this->hap2 << "\t" << this->start_site << "\t" << this->end_site << "\t" << this->len_cm << std::endl;
        }

        bool operator==(const Match &other){
            if((this->hap1 == other.hap1) && (this->hap2 == other.hap2) && (this->start_site == other.start_site) && (this->end_site == other.end_site)){
                return true;
            }
            else{
                return false;
            }
        }
};

int nextStart(int hap1, int hap2, int start, int max_gap, std::vector<int> &site_mapping, std::vector<Match> &matches, rateMapData &gen_map, float min_seed, float min_extend, int min_seed_markers, int min_extend_markers, std::vector<std::vector<int>> &genotype_array);
int extendStart(int hap1, int hap2, int start, int max_gap, std::vector<int> &site_mapping, std::vector<Match> &matches, rateMapData &gen_map, float min_seed, float min_extend, int min_seed_markers, int min_extend_markers, std::vector<std::vector<int>> &genotype_array);
int nextInclEnd(int hap1, int hap2, int incl_end, int max_gap, std::vector<int> &site_mapping, std::vector<Match> &matches, rateMapData &gen_map, float min_seed, float min_extend, int min_extend_markers, std::vector<std::vector<int>> &genotype_array);
int extendInclEnd(int hap1, int hap2, int incl_end, int max_gap, std::vector<int> &site_mapping, std::vector<Match> &matches, rateMapData &gen_map, float min_seed, float min_extend, int min_extend_markers, std::vector<std::vector<int>> &genotype_array);
std::string hapToTskId(int hap);
std::string processSeed(int hap1, int hap2, int start, int incl_end, int max_gap, std::vector<int> &site_mapping, std::vector<Match> &matches, rateMapData &gen_map, float min_seed, float min_extend, int min_seed_markers, int min_extend_markers, float min_output, std::vector<std::vector<int>> &genotype_array);



