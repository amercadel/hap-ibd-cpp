## hap-ibd-cpp


#### Intro
This is an implemntation of the Browning Lab's hap-ibd, an algorithm for the detection of segments that are identical by descent. As of right now, it is not optimized and still quite slow, but as of 10/28/24, the output will be same, with the exception of the homozygous by descent segments and identical by descent segments being included in the same file.

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
