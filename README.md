## hap-ibd-cpp


#### Intro
This is an implemntation of the Browning Lab's hap-ibd, an algorithm for the detection of segments that are identical by descent. The implementation is almost exactly the same with a marginal slowdown. This is written in C++ for ease of integration of other tools that use IBD analysis. 

Original paper can be found here: 
Zhou, Y., Browning, S. R., & Browning, B. L. (2020). A Fast and Simple Method for Detecting Identity-by-Descent Segments in Large-Scale Data. American journal of human genetics, 106(4), 426–437. https://doi.org/10.1016/j.ajhg.2020.02.010

#### Installation
##### Requirements
- meson (1.5.2)
- ninja (1.12.1)

##### Setup
- After cloning the repo, cd into hap-ibd-cpp, then run: ```git submodule update --init --recursive```
- Once you've done that, return to the root directory and run: ```meson setup build```
- Then, cd into the build directory, and run: ```ninja```
- This will create the hap-ibd-cpp executable, which for now only takes in two arguments: a vcf file and a plink genetic map (need to fix to accept nondefult parameters)
- The vcf file must not have any multiallelic sites (can be done using bcftools)
- ex: ```./hap-ibd-cpp {vcf file path} {plink rate map file path}```
- The output will be a file called output.txt


#### Future
- may look into exposing the bindings to python, but this would require a big cleanup of the code