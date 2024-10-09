#pragma once
#include <vector>
#include <iostream>
#include <string>
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


std::vector<std::string> splitVCFByPos(const char* input_vcf, std::vector<int> overlapping_windows){
    std::vector<std::string> split_files;
    for(size_t i = 0; i < overlapping_windows.size(); i = i + 2){
        int curr_index = 0;
        htsFile *input = hts_open(input_vcf, "r");
        bcf_hdr_t *hdr = bcf_hdr_read(input);
        bcf1_t *rec = bcf_init();
        int start_pos = overlapping_windows[i];
        int end_pos = overlapping_windows[i + 1];
        char output_filename[256];
        sprintf(output_filename, "intermediate_vcf_%d_%d.vcf", start_pos, end_pos);
        split_files.push_back(std::string(output_filename));
        htsFile* out = hts_open(output_filename, "w");
        bcf_hdr_write(out, hdr);
        
        while(bcf_read(input, hdr, rec) == 0){
            if (curr_index < start_pos){
                curr_index++;
            }
            else if(curr_index > end_pos){
                hts_close(out);
                hts_close(input);
                bcf_destroy(rec);
                bcf_hdr_destroy(hdr);
                break;
            }
            else{
                bcf_write(out, hdr, rec);
                curr_index++;
            }
        }
    }
    

    

}