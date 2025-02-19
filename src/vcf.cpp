#include <stdexcept>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
extern "C"
{
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/vcfutils.h"
}
#include "vcf.hpp"
#include "utils.hpp"



void getSiteMappingAndGenotypes(char* vcf_file, std::vector<std::vector<int>> &genotype_array, std::vector<int> &site_mapping, int n_threads){
    
    htsFile *fp = hts_open(vcf_file, "r");
    if(!fp){
        std::cerr << "Error: htslib is unable to open the VCF file\n";
        hts_close(fp);
    }
    bcf_hdr_t *hdr = bcf_hdr_init("r");
    hdr = bcf_hdr_read(fp);
    if(!hdr){
        std::cerr << "Error, htslib failed to read VCF Header";
        bcf_hdr_destroy(hdr);
        hts_close(fp);
    }
    bcf1_t *rec = bcf_init();
    if (!rec) {
        std::cerr << "Error allocating memory for VCF record.\n";
        bcf_destroy(rec);
        bcf_hdr_destroy(hdr);
        hts_close(fp);
    }

    int n_samples = bcf_hdr_nsamples(hdr);
    

   
    int *gt_arr = NULL, ngt_arr = 0;
    while(bcf_read(fp, hdr, rec) == 0){
        int32_t pos = rec->pos + 1;
        site_mapping.push_back(pos);
        bcf_unpack(rec, BCF_UN_ALL);
        std::vector<int> alleles;
        alleles.resize(2 * n_samples);



        if(bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr) > 0){
            for(int i = 0; i < n_samples; i++){
                int allele1 = gt_arr[2*i];
                int allele1_gt = bcf_gt_allele(allele1);
                int allele2 = gt_arr[2*i + 1];
                int allele2_gt = bcf_gt_allele(allele2);
                alleles[2 * i] = allele1_gt;
                alleles[2 * i + 1] = allele2_gt;
                
            }
        }
        genotype_array.emplace_back(std::move(alleles));
        
        
    }
    free(gt_arr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
}


