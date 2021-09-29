hmcnc_cpp_sources = files([
                  'SamToBed.cpp',
                  'viterbi.cpp',
                  ])

 executable(
     'hmcnc',
     files([
        'hmcnc_defs.cpp',
        'hmcnc_io.cpp',
        'hmmcnc.cpp',
        'main.cpp'
         ]),
     install : true,
     dependencies : hmcnc_lib_deps,
     cpp_args : hmcnc_cpp_args,
     )
