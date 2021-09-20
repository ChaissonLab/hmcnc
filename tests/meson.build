hmmcnc_tests_cpp = files([
  'test_hmmcnc.cpp',
  '../src/hmmcnc.cpp'
])

hmmcnc_gtest_dep = dependency(
  'gtest',
  main: true,
  fallback: ['gtest', 'gtest_main_dep']
)

hmmcnc_test = executable(
  'hmmcnc_test',
  hmmcnc_tests_cpp,
  dependencies : [hmcnc_lib_deps, hmmcnc_gtest_dep],
  cpp_args : hmcnc_cpp_args,
  install : false
)

test(
  'hmmcnc unit tests',
  hmmcnc_test,
)
