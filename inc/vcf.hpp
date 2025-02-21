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

// populates the genotype array and site mapping vectors 
void getSiteMappingAndGenotypes(char* vcf_file, std::vector<std::vector<int>> &genotype_array, std::vector<int> &site_mapping);
void getSiteMappingAndGenotypes(char* vcf_file, std::vector<std::unordered_set<int>> &alternate_allele_map, std::vector<int> &site_mapping);

// helper function for the hashset implementation. Needed to fetch genotype information since it is not all stored in memory
// OUTPUT: int (0 or 1) representing reference allele or alternate allele
int getHaplotype(int site_index, int haplotype_index, std::vector<std::unordered_set<int>> &alt_map);

#endif