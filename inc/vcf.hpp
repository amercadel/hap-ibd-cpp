#pragma once
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <bitset>


#define MAX_N_SAMPLES 5000







void getSiteMappingAndGenotypes(char* vcf_file, std::vector<std::vector<int>> &genotype_array, std::vector<int> &site_mapping, int n_threads);
void getSiteMappingAndGenotypes(char* vcf_file, std::unordered_map<int, std::bitset<MAX_N_SAMPLES>> &alt_map, std::vector<int> &site_mapping, int n_threads);
std::vector<char*> splitVCFByPos(char* input_vcf, std::vector<std::pair<int, int>> &overlapping_windows);
int getAllele(int site, int hap, std::unordered_map<int, std::bitset<MAX_N_SAMPLES>> &alt_map);
