project(
        'hmcnc',
        ['cpp', 'c'],
        version : '0.0.1',
        default_options : [
                        'buildtype=release',
                        'warning_level=3',
                        'cpp_std=c++14',
                        'c_std=c11',
                        'b_ndebug=if-release'],
       license : 'MIT',
       meson_version : '>= 0.52.0',
       )

############
# CXXFLAGS #
############

hmcnc_warning_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
        '-Wduplicated-cond',
        '-Wduplicated-branches',
        '-Wlogical-op',
        '-Wrestrict',
        '-Wnull-dereference',
        '-Wuseless-cast',
        '-Wdouble-promotion',
        '-Wshadow',
        '-Wformat=1',
        ]
        if cpp.has_argument(cflag)
           hmcnc_warning_flags += cflag
        endif
endforeach

hmcnc_cpp_args = hmcnc_warning_flags

################
# dependencies #
################

# boost
boost_dep = dependency('boost', required : true, include_type : 'system')

# threads
threads_dep = dependency('threads', required : true)

# htslib
htslib_dep = dependency('htslib', required : true, version : '>=1.4', fallback : ['htslib', 'htslib_dep'])

# zlib
zlib_dep = dependency('zlib', required : true, fallback : ['zlib', 'zlib_dep'])

hmcnc_lib_deps = [
               boost_dep,
               threads_dep,
               zlib_dep,
               htslib_dep,
               ]

#####################
# sources + library #
#####################

subdir('src')
subdir('tests')
