#pragma once
#include <vector>
#include <iostream>
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
