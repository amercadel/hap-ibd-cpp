#ifndef VCF_HPP
#define VCF_HPP
#include <vector>
#include <unordered_set>
#include <stdexcept>
#include <fstream>

extern "C"
{
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/vcfutils.h"
}


void getSiteMappingAndGenotypes(char* vcf_file, std::vector<std::vector<int>> &genotype_array, std::vector<int> &site_mapping);
void getSiteMappingAndGenotypes(char* vcf_file, std::vector<std::unordered_set<int>> &alternate_allele_map, std::vector<int> &site_mapping);


int getHaplotype(int site_index, int haplotype_index, std::vector<std::unordered_set<int>> &alt_map);

#endif