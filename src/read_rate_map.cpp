#include <sstream>
#include <fstream>
#include <string>
#include "read_rate_map.hpp"
#include "utils.hpp"

double getGeneticPosition(std::vector<double> &interpolated_cm, int site_idx){
    // assert(site_idx < interpolated_cm.size());
    // if(site_idx < interpolated_cm.size()){
    //     return interpolated_cm[site_idx];
    // }
    // else{
    //     std::cerr << "Site outside of interolated map\n";
    //     exit(EXIT_FAILURE);
    // }
    return interpolated_cm[site_idx];
	
}

std::vector<std::string> split(std::string &line, char delim){
    std::vector<std::string> split_string;
    std::string str;
    std::stringstream ss(line);
    while (std::getline(ss, str, delim)){
        split_string.push_back(str);
    }
    return split_string;

}

std::vector<double> rateMapData::interpolateVector(std::vector<int> &sites) {
    std::vector<double> interpolated_cm;
    for(size_t c = 0; c < sites.size(); c++){
        interpolated_cm.push_back(genPos(sites[c]));
    }
    return interpolated_cm;
}

double rateMapData::genPos(int site){
    int map_size_m1 = bp_vec.size() - 1;
    double min_end_cm_dist = 5.0f;
    auto it = std::lower_bound(bp_vec.begin(), bp_vec.end(), site);
    int index = it - bp_vec.begin();
    if (it != bp_vec.end() && *it == site){
        return cm_vec[index];
    }else{
        int insertion_point = index;
        int a_index = insertion_point - 1;
        int b_index = insertion_point;
        if (a_index == map_size_m1){
            auto gp_it = std::lower_bound(cm_vec.begin(), cm_vec.end(), cm_vec.back() - min_end_cm_dist);
            insertion_point = gp_it - cm_vec.begin();
            if(insertion_point < 0){
                insertion_point = -insertion_point - 2;
            }
            a_index = std::max(insertion_point, 0);
            b_index = map_size_m1;
        }
        else if(b_index == 0){
            auto gp_it = std::lower_bound(cm_vec.begin(), cm_vec.end(), cm_vec.front() + min_end_cm_dist);
            insertion_point = gp_it - cm_vec.begin();
            if(insertion_point < 0){
                insertion_point = -insertion_point - 2;
            }
            a_index = 0;
            b_index = std::min(insertion_point, map_size_m1);
        }
        int x = site;
        int a = bp_vec[a_index];
        int b = bp_vec[b_index];
        double fa = cm_vec[a_index];
        double fb = cm_vec[b_index];
        return fa + (((double)(x - a) / (double)(b - a)) * (fb - fa));
    }
    
}

double rateMapData::interpolateBasePairToGenPos(int site){
    int idx = findVectorIndex(bp_vec, site);
    if (idx > 0){
        return cm_vec[idx];
    }
    else{
        idx = findInsertionIndex(bp_vec, site);
        if (idx == 0){
            return cm_vec[0];
        }
        else if(idx >= bp_vec.size()){
            double slope = ((cm_vec.back() - cm_vec.front()) / (bp_vec.back() - bp_vec.front()));
            double y_int = cm_vec.back() - (slope * bp_vec.back());
            return slope * site + y_int;
        }

        else{
            int x0 = bp_vec[idx - 1];
            int x1 = bp_vec[idx];
            double y0 = cm_vec[idx - 1];
            double y1 = cm_vec[idx];
            double y = ((y0 * (x1 - site)) + (y1 * (site - x0))) / (x1 - x0);
            return y;

        }
    }
}


rateMapData readRateMap(char* filename, std::vector<int> &sites){
    std::ifstream inputFile;
    std::string line;
    std::vector<int> bp_vec;
    std::vector<double> cm_vec;
    inputFile.open(filename);
    rateMapData res;
    while (getline(inputFile, line)) {
        int bp;
        double cm;
        std::string tmp;
        std::stringstream input_str(line); 
        std::vector<std::string> data = split(line, ' ');
        bp = std::stoi(data[3]);
        cm = std::stod(data[2]);
        bp_vec.push_back(bp);
        cm_vec.push_back(cm);

    }
    
    res.bp_vec = bp_vec;
    res.cm_vec = cm_vec;
    res.last_bp = bp_vec[bp_vec.size() - 1];
    res.last_cm = cm_vec[cm_vec.size() - 1];
    res.interpolated_cm = res.interpolateVector(sites);
    return res;
}
