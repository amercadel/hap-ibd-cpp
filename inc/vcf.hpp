#pragma once
#include <vector>






void getSiteMappingAndGenotypes(char* vcf_file, std::vector<std::vector<int>> &genotype_array, std::vector<int> &site_mapping);
std::vector<char*> splitVCFByPos(char* input_vcf, std::vector<std::pair<int, int>> &overlapping_windows);