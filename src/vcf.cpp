#include "vcf.hpp"
#include "utils.hpp"



void getSiteMappingAndGenotypes(char* vcf_file, std::vector<std::vector<int>> &genotype_array, std::vector<int> &site_mapping){
    
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




void getSiteMappingAndGenotypes(char* vcf_file, std::vector<std::unordered_set<int>> &alternate_allele_map, std::vector<int> &site_mapping){
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
        site_mapping.emplace_back(pos);
        bcf_unpack(rec, BCF_UN_ALL);
        std::unordered_set<int> alternate_alleles;
        alternate_alleles.reserve(n_samples * 2);

        if(bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr) > 0){
            for(int i = 0; i < n_samples; i++){
                int allele1 = gt_arr[2*i];
                int allele1_gt = bcf_gt_allele(allele1);
                int allele2 = gt_arr[2*i + 1];
                int allele2_gt = bcf_gt_allele(allele2);
                if (allele1_gt == 1){
                    alternate_alleles.insert(2 * i);
                }
                if (allele2_gt == 1){
                    alternate_alleles.insert(2 * i + 1);
                }
            }
        }
        alternate_allele_map.emplace_back(std::move(alternate_alleles));
        
        
    }
    free(gt_arr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
}



int getHaplotype(int site_index, int haplotype_index, std::vector<std::unordered_set<int>> &alt_map){
    std::unordered_set<int> &site = alt_map[site_index];
    if(site.find(haplotype_index) != site.end()){
        return 1;
    } else {
        return 0; 
    }
}
