#
# GTest
#
hmcnc_gtest_dep = dependency(
  'gtest',
  main: true,
  fallback: ['gtest', 'gtest_main_dep']
)

#
# Test runner
#
hmcnc_tests_cpp_sources = hmcnc_cpp_sources
hmcnc_tests_cpp_sources += files([
  'test_hmcnc.cpp',
  'test_hmcnc_io.cpp',
])

hmcnc_test = executable(
  'hmcnc_test',
  hmcnc_tests_cpp_sources,
  dependencies : [hmcnc_lib_deps, hmcnc_gtest_dep],
  cpp_args : hmcnc_cpp_args,
  install : false
)

test(
  'hmcnc unit tests',
  hmcnc_test,
)
