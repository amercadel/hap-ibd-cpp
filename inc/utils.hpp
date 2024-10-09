#include <vector>


template<typename T>
T roundToNDigits(T num, int n_digits){
    auto scale = std::pow(10.0, n_digits);
    return std::round(num * scale) / scale;
}

template<typename T>
int findVectorIndex(std::vector<T> &vec, T val){
    int left = 0;
    int right = vec.size() - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (vec[mid] == val) {
            return mid;
        } else if (vec[mid] < val) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    return -1;
}
template<typename T>
int findInsertionIndex(std::vector<T> &vec, T val){
    int left = 0;
    int right = vec.size() - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (vec[mid] == val) {
            return mid;
        } else if (vec[mid] < val) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    return left;
}
template<typename T>
int findInsertionIndex(std::vector<T> &vec, int from, int to, T val){
    int left = from;
    int right = to;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (vec[mid] == val) {
            return mid;
        } else if (vec[mid] < val) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    return left;
}


std::vector<int> overlappingWindows(std::vector<float> cm, float min_seed, int min_markers, int n_threads){
    std::vector<int> starts;
    std::vector<int> ends;
    float L = cm.back() - cm.front();
    float window_size = ((L - min_seed) / n_threads + min_seed);
    float dist = (L - min_seed) / n_threads;
    float start_gen_pos = 0.0f;
    int i = 0;
    while(i < n_threads){
        float end_gen_pos = start_gen_pos + window_size;
        starts.push_back(findInsertionIndex(cm, start_gen_pos));
        ends.push_back(findInsertionIndex(cm, end_gen_pos));
        start_gen_pos += dist;
        i++;
    }
    std::vector<int> window_vec;
    for(size_t i = 0; i < starts.size(); i++){
        window_vec.push_back(starts[i]);
        window_vec.push_back(ends[i]);
        std::cout << starts[i] << " " << ends[i] << std::endl;
    }
    return window_vec;
}


