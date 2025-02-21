project('test-hap-ibd-cpp', ['c', 'cpp'], default_options: ['c_std=c99', 'cpp_std=c++20', 'warning_level=0'])

zdep = dependency('zlib', version : '>=1.2.8')
htslib_dep = dependency('htslib', required: true)

pbwt_sources = files([
    'pbwt_am/array.c',
    'pbwt_am/hash.c',
    'pbwt_am/dict.c',
    'pbwt_am/pbwtCore.c',
    'pbwt_am/pbwtIO.c',
    'pbwt_am/pbwtMatch.c',
    'pbwt_am/pbwtSample.c',
    'pbwt_am/utils.c',
    'pbwt_am/pbwtHtslib.c',
    'pbwt_am/array.h',
    'pbwt_am/hash.h',
    'pbwt_am/dict.h',
    'pbwt_am/utils.h', 
    'pbwt_am/pbwtImpute.c'
])

pbwt_inc = include_directories('pbwt_am')
pbwt_dep = [htslib_dep]

pbwt_lib = static_library('pbwt', pbwt_sources, dependencies:  [pbwt_dep, zdep], include_directories: pbwt_inc)

pbwt_dep = declare_dependency(include_directories: pbwt_inc, 
                            link_with: pbwt_lib)


my_includes = include_directories('inc')

executable('hap-ibd-cpp', ['src/main.cpp', 'src/match.cpp', 'src/vcf.cpp', 'src/read_rate_map.cpp', 'src/utils.cpp'], dependencies: [htslib_dep, pbwt_dep], include_directories: [my_includes], override_options : ['c_std=c20'])
