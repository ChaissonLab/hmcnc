##############
# CRAM tests #
##############

hmcnc_cram_script = find_program('cram')
hmcnc_cram_tests = [
  'chr22',
]

foreach t : hmcnc_cram_tests
  test(
    'hmcnc cram test - ' + t,
    hmcnc_cram_script,
    args : [
      '--xunit-file=' + join_paths(meson.build_root(), 'hmcnc-cram-' + t.underscorify() + '-test.xml'),
      '--verbose'] + files(t + '.t'),
    env : [
      '__HMCNC_EXE=' + hmcnc_exe.full_path(),
    ],
    timeout : 86400
  )
endforeach
