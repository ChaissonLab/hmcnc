hmcnc_gtest_dep = dependency(
  'gtest',
  main: true,
  fallback: ['gtest', 'gtest_main_dep']
)

hmcnc_unit_test_sources = files([
  'test_hmcnc.cpp',
  'test_hmcnc_io.cpp',
])
hmcnc_unit_test_sources += hmcnc_cpp_sources

hmcnc_test = executable(
  'hmcnc_test',
  hmcnc_unit_test_sources,
  dependencies : [hmcnc_lib_deps, hmcnc_gtest_dep],
  cpp_args : hmcnc_cpp_args,
  include_directories: hmcnc_include_directories,
  install : false
)

test(
  'hmcnc unit tests',
  hmcnc_test,
)
