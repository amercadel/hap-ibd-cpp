#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <unistd.h>
#include "match.hpp"
#include "read_rate_map.hpp"
#include "vcf.hpp"
extern "C"
{
#include "pbwt.h"
}

int main()