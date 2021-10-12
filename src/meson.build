hmcnc_cpp_sources = files([
  'SamToBed.cpp',
  'viterbi.cpp',
])

hmcnc_cpp_sources = files([
  'hmmcnc.cpp',
  'hmcnc_io.cpp',
])

hmcnc_exe_sources = files([
  'main.cpp',
])
hmcnc_exe_sources += hmcnc_cpp_sources

hmcnc_exe = executable(
  'hmcnc',
  hmcnc_exe_sources,
  install : true,
  dependencies : hmcnc_lib_deps,
  cpp_args : hmcnc_cpp_args,
)
