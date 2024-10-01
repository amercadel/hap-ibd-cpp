## hap-ibd-cpp

This is an implemntation of the Browning Lab's hap-ibd, an algorithm for the detection of segments that are identical by descent.

I use meson and ninja as my build system, so make sure you have those installed.

After cloning the repo, cd into hap-ibd-cpp, then run git submodule update --init --recursive

Once you've done that, run meson setup build.

Then, cd into the build directory, and run ninja

This will create the hap-ibd-cpp executable, which for now only takes in two arguments: a vcf file and a plink genetic map

The vcf file must not have any multiallelic sites

ex: ./hap-ibd-cpp {vcf file path} {plink rate map file path}

The output will be a file called output_matches.txt
