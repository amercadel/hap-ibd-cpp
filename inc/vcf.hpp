#pragma once
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include "utils.hpp"

extern "C"
{
#include "htslib/vcf.h"
#include "htslib/hts.h"
}




std::vector<int> getSiteMapping(char* vcf_file){
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
    std::vector<int> ret;
    while(bcf_read(fp, hdr, rec) == 0){
        int32_t pos = rec->pos + 1;
        ret.push_back(pos);
    }
    
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    return ret;
}


std::vector<char*> splitVCFByPos(char* input_vcf, std::vector<std::pair<int, int>> overlapping_windows){
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
        std::strcpy(output_file_name_cstr, output_file_name.c_str());
        file_names.push_back(output_file_name_cstr);
        output_files[i] = hts_open(output_file_name_cstr, "w");  // Open output file
        bcf_hdr_write(output_files[i], hdr);
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