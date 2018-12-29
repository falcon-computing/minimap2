FILE(REMOVE_RECURSE
  "CMakeFiles/googletest-download"
  "CMakeFiles/googletest-download-complete"
  "deps/googletest/src/googletest-download-stamp/googletest-download-install"
  "deps/googletest/src/googletest-download-stamp/googletest-download-mkdir"
  "deps/googletest/src/googletest-download-stamp/googletest-download-download"
  "deps/googletest/src/googletest-download-stamp/googletest-download-update"
  "deps/googletest/src/googletest-download-stamp/googletest-download-patch"
  "deps/googletest/src/googletest-download-stamp/googletest-download-configure"
  "deps/googletest/src/googletest-download-stamp/googletest-download-build"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/googletest-download.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
