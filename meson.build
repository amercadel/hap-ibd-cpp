project('test-hap-ibd-cpp', ['c', 'cpp'], default_options: ['c_std=c99', 'cpp_std=c++17', 'warning_level=0'])

zdep = dependency('zlib', version : '>=1.2.8')
htslib_dep = dependency('htslib', required: true)

pbwt_sources = files([
    'pbwt/array.c',
    'pbwt/hash.c',
    'pbwt/dict.c',
    'pbwt/pbwtCore.c',
    'pbwt/pbwtIO.c',
    'pbwt/pbwtMatch.c',
    'pbwt/pbwtSample.c',
    'pbwt/utils.c',
    'pbwt/pbwtHtslib.c',
    'pbwt/array.h',
    'pbwt/hash.h',
    'pbwt/dict.h',
    'pbwt/utils.h', 
    'pbwt/pbwtImpute.c'
])

pbwt_inc = include_directories('pbwt')
pbwt_dep = [htslib_dep]

pbwt_lib = static_library('pbwt', pbwt_sources, dependencies:  [pbwt_dep, zdep], include_directories: pbwt_inc)

pbwt_dep = declare_dependency(include_directories: pbwt_inc, 
                            link_with: pbwt_lib)


my_includes = include_directories('inc')

executable('hap-ibd-cpp', ['src/main.cpp', 'src/match.cpp', 'src/vcf.cpp', 'src/read_rate_map.cpp', 'src/utils.cpp'], dependencies: [htslib_dep, pbwt_dep], include_directories: [my_includes])
