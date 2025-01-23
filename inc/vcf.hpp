#pragma once
#include <vector>






void getSiteMappingAndGenotypes(char* vcf_file, std::vector<std::vector<int>> &genotype_array, std::vector<int> &site_mapping, int n_threads);
void getSiteMappingAndGenotypes(char* vcf_file, std::map<int, std::vector<int>> &alt_map, std::vector<int> &site_mapping, int n_threads);
std::vector<char*> splitVCFByPos(char* input_vcf, std::vector<std::pair<int, int>> &overlapping_windows);
int getAllele(int site, int hap, std::map<int, std::vector<int>> &alt_map);
// void new_genotype_array(char* vcf_file);