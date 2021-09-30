hmcnc_tests_cpp = files([
  'test_hmcnc.cpp',
  '../src/hmmcnc.cpp',
  '../src/hmcnc_io.cpp',
])

hmcnc_gtest_dep = dependency(
  'gtest',
  main: true,
  fallback: ['gtest', 'gtest_main_dep']
)

hmcnc_test = executable(
  'hmcnc_test',
  hmcnc_tests_cpp,
  dependencies : [hmcnc_lib_deps, hmcnc_gtest_dep],
  cpp_args : hmcnc_cpp_args,
  install : false
)

test(
  'hmcnc unit tests',
  hmcnc_test,
)
