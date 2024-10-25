#pragma once
#include <iostream>
#include <vector>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <string>
#include "read_rate_map.hpp"
 




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
static bool compareMatch(const Match &m1, const Match &m2) {
    if (m1.hap1 != m2.hap1) {
        return m1.hap1 < m2.hap1;
    }
    if (m1.hap2 != m2.hap2) {
        return m1.hap2 < m2.hap2;
    }
    if (m1.start_site != m2.start_site) {
        return m1.start_site < m2.start_site;
    }
    return m1.end_site < m2.end_site;
}

// static bool compareStartSite(const Match &m1, const Match &m2){
//     // m1 must be to the left of m2
//     return m1.start_site < m2.start_site;
// }

// static bool compareHaplotype1(const Match &m1, const Match &m2){
//     return m1.hap1 < m2.hap1;
// }

// static bool compareHaplotype2(const Match &m1, const Match &m2){
//     return m1.hap2 < m2.hap2;
// }

// static bool compareHaps(const Match &m1, const Match &m2) {
//     if (m1.hap1 == m2.hap1) {
//         return m1.hap2 < m2.hap2;
//     }
//     return m1.hap1 < m2.hap1;
// }


// bool checkAllele(int hap1, int hap2, int site, std::vector<Match> &matches){
//     // matches should be sorted by hap1 and hap2
//     int idx = 0;
//     while ((matches[idx].hap1 != hap1) && (matches[idx].hap2 != hap2) && (idx < matches.size())){
//         idx++;
//     }
//     while((matches[idx].hap1 == hap1) && (matches[idx].hap2 == hap2) && (idx < matches.size())){
//         if((matches[idx].start_site <= site) && (site <= matches[idx].end_site)){
//             return true;
//         }
//         idx++;
//     }
//     return false;
// }


int nextStart(int hap1, int hap2, int start, int max_gap, std::vector<int> &site_mapping, std::vector<Match> &matches, rateMapData &gen_map, float min_seed, float min_extend, int min_seed_markers, int min_extend_markers, std::vector<std::vector<int>> &genotype_array){
    if (start < 2 || max_gap < 0){
        return start;
    }
    int m = start - 1;
    int first_mismatch_pos = site_mapping[m];

    int first_match = start - 2;
    while(m > 0){
        --m;
        int a1 = genotype_array[m][hap1];
        int a2 = genotype_array[m][hap2];
        // std::cout << m << " " << a1 << " " << a2 << std::endl;
        if(a1!=a2){
            
            if((first_mismatch_pos - site_mapping[m]) > max_gap){
                ++m;
                break;
            }
            else if (m > 0){
                first_match = m - 1;
            }
        }
    }
    
    float len = gen_map.interpolated_cm[first_match] - gen_map.interpolated_cm[m];
    if (len >= min_seed && ((first_match - m) >= min_seed_markers - 1)){
        return -1;
    }
    else{
        int ret = (len < min_extend || ((first_match - m) < min_extend_markers - 1)) ? start : m;
        return ret;
    }
}


int extendStart(int hap1, int hap2, int start, int max_gap, std::vector<int> &site_mapping, std::vector<Match> &matches, rateMapData &gen_map, float min_seed, float min_extend, int min_seed_markers, int min_extend_markers, std::vector<std::vector<int>> &genotype_array){
    int prev_start = start;
    int next_start = nextStart(hap1, hap2, prev_start, max_gap, site_mapping, matches, gen_map, min_seed, min_extend, min_seed_markers, min_extend_markers, genotype_array);
    while (next_start>=0 && (next_start < prev_start)) {
        prev_start = next_start;
        next_start = nextStart(hap1, hap2, prev_start, max_gap, site_mapping, matches, gen_map, min_seed, min_extend, min_seed_markers, min_extend_markers, genotype_array);
    }
    return next_start;
}

int nextInclEnd(int hap1, int hap2, int incl_end, int max_gap, std::vector<int> &site_mapping, std::vector<Match> &matches, rateMapData &gen_map, float min_seed, float min_extend, int min_extend_markers, std::vector<std::vector<int>> &genotype_array){
    int last_marker = site_mapping.size() - 1;
    if ((incl_end>(last_marker - 2)) || (max_gap < 0)) {
        return incl_end;
    }
    int m = incl_end + 1;
    int first_mismatch_pos = site_mapping[m];
    int first_match = incl_end + 2;
    while(m < last_marker){
        ++m;
        int a1 = genotype_array[m][hap1];
        int a2 = genotype_array[m][hap2];
        if(a1 != a2){
            if ((site_mapping[m] - first_mismatch_pos) > max_gap){
                --m;
                break;
            }
            else if(m < last_marker){
                first_match = m + 1;
            }
        }

    }
    float len = (gen_map.interpolated_cm[m] - gen_map.interpolated_cm[first_match]);
    return (len<min_extend || (m-first_match)<(min_extend_markers - 1)) ? incl_end : m;
}

int extendInclEnd(int hap1, int hap2, int incl_end, int max_gap, std::vector<int> &site_mapping, std::vector<Match> &matches, rateMapData &gen_map, float min_seed, float min_extend, int min_extend_markers, std::vector<std::vector<int>> &genotype_array){
    int last_marker = site_mapping.size();
    while (incl_end<last_marker && (genotype_array[incl_end + 1][hap1] == genotype_array[incl_end + 1][hap2])) {
        ++incl_end;
    }
    int prev_incl_end = incl_end;
    int next_incl_end = nextInclEnd(hap1, hap2, prev_incl_end, 1000, site_mapping, matches, gen_map, 2.0, 1.0, 100, genotype_array);
    while (next_incl_end>prev_incl_end) {
        prev_incl_end = next_incl_end;
        next_incl_end = nextInclEnd(hap1, hap2, prev_incl_end, 1000, site_mapping, matches, gen_map, 2.0, 1.0, 100, genotype_array);
    }
    return next_incl_end;
}


std::string hapToTskId(int hap){
    int id = hap / 2;
    int hap_id = (hap % 2) + 1;
    std::string ret;
    ret = "tsk_" + std::to_string(id) + "\t" + std::to_string(hap_id);
    return ret;
}


std::string processSeed(int hap1, int hap2, int start, int incl_end, int max_gap, std::vector<int> &site_mapping, std::vector<Match> &matches, rateMapData &gen_map, float min_seed, float min_extend, int min_seed_markers, int min_extend_markers, float min_output, std::vector<std::vector<int>> &genotype_array){
    std::stringstream out;
    
    start = extendStart(hap1, hap2, start, max_gap, site_mapping, matches, gen_map, min_seed, min_extend, min_seed_markers, min_extend_markers, genotype_array);
    if (start>=0) {
        incl_end = extendInclEnd(hap1, hap2, incl_end, max_gap, site_mapping, matches, gen_map, min_seed, min_extend, min_extend_markers, genotype_array);
        if ((gen_map.interpolated_cm[incl_end] - gen_map.interpolated_cm[start])>=min_output) {
            out << hapToTskId(hap1) << "\t" << hapToTskId(hap2) << "\t" << "20" << "\t" << site_mapping[start] << "\t" << site_mapping[incl_end] << "\t" << roundToNDigits(gen_map.interpolated_cm[incl_end] - gen_map.interpolated_cm[start], 3) << "\n";
            }
        return out.str();
        }
    else{
        return "";
    }
}



