

struct hap_ibd_parameters {
    char* input_vcf;
    char* map_file;
    char* output_file;
    double min_seed = 2.0;
    int max_gap = 1000;
    double min_extend = 1.0;
    double min_output = 2.0;
    int min_mac = 2;
    int min_markers = 100;
    int n_threads = 1;
};