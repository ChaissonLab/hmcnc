hmcnc_cpp_sources = files([
                  'SamToBed.cpp',
                  'viterbi.cpp',
                  ])

 executable(
     'hmcnc',
     files([
        'hmmcnc.cpp',
         ]),
     install : true,
     dependencies : hmcnc_lib_deps,
     cpp_args : hmcnc_cpp_args,
     )