#include <stdexcept>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
extern "C"
{
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/vcfutils.h"
}
#include "vcf.hpp"
#include "utils.hpp"



void getSiteMappingAndGenotypes(char* vcf_file, std::vector<std::vector<int>> &genotype_array, std::vector<int> &site_mapping){
    htsFile *fp = hts_open(vcf_file, "r");
    if(!fp){
        std::cerr << "Error: htslib is unable to open the VCF file\n";
    }
    bcf_hdr_t *hdr = bcf_hdr_init("r");
    hdr = bcf_hdr_read(fp);
    bcf1_t *rec = bcf_init();
    if (!rec) {
        std::cerr << "Error allocating memory for VCF record.\n";
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
    }
    while(bcf_read(fp, hdr, rec) == 0){
        int32_t pos = rec->pos + 1;
        site_mapping.push_back(pos);
        bcf_unpack(rec, BCF_UN_ALL);
        std::vector<int> alleles;

        int n_samples = bcf_hdr_nsamples(hdr);
        int *gt_arr = NULL, ngt_arr = 0;
        if(bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr) > 0){
            for(int i = 0; i < n_samples; i++){
                int allele1 = gt_arr[2*i];
                int allele1_gt = bcf_gt_allele(allele1);
                int allele2 = gt_arr[2*i + 1];
                int allele2_gt = bcf_gt_allele(allele2);
                alleles.push_back(allele1_gt);
                alleles.push_back(allele2_gt);
                
            }
        }
        genotype_array.push_back(alleles);
        free(gt_arr);
        
    }
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
}




std::vector<char*> splitVCFByPos(char* input_vcf, std::vector<std::pair<int, int>> &overlapping_windows){
    htsFile *input = hts_open(input_vcf, "r");
    bcf_hdr_t *hdr = bcf_hdr_read(input);
    bcf1_t *rec = bcf_init();
    int curr_index = 0;
    std::vector<char*> file_names;
    std::vector<htsFile*> output_files(overlapping_windows.size());
    for(size_t i = 0; i < output_files.size(); i++){
        std::ostringstream oss;
        oss << "intermediate_vcf_" << overlapping_windows[i].first << "_" << overlapping_windows[i].second << ".vcf";
        std::string output_file_name = oss.str();
        char* output_file_name_cstr = new char[output_file_name.length() + 1];
        strcpy(output_file_name_cstr, output_file_name.c_str());
        file_names.push_back(output_file_name_cstr);
        output_files[i] = hts_open(output_file_name_cstr, "w");  // Open output file
        int ret = bcf_hdr_write(output_files[i], hdr);
        if (ret != 0){
            std::cerr << "failed to write new vcf file\n";
            exit(EXIT_FAILURE);
        }

    }
    while(bcf_read(input, hdr, rec) == 0){
        for(size_t c = 0; c < output_files.size(); c++){ // iterate through output files
            if(overlapping_windows[c].first <= curr_index && curr_index <= overlapping_windows[c].second){ // check if current index is between the corresponding window
                int ret = bcf_write(output_files[c], hdr, rec); // if so, write that line/record to the current file
                if(ret == 1){
                    std::cerr << "HTSLIB failed to write record";
                }
            }
        }
        curr_index++;
    }
    for(size_t i = 0; i < output_files.size(); i++){
        hts_close(output_files[i]);
    }
    
    hts_close(input);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    return file_names;
    

}