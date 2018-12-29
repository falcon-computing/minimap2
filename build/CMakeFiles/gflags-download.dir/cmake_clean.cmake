FILE(REMOVE_RECURSE
  "CMakeFiles/gflags-download"
  "CMakeFiles/gflags-download-complete"
  "deps/gflags/src/gflags-download-stamp/gflags-download-install"
  "deps/gflags/src/gflags-download-stamp/gflags-download-mkdir"
  "deps/gflags/src/gflags-download-stamp/gflags-download-download"
  "deps/gflags/src/gflags-download-stamp/gflags-download-update"
  "deps/gflags/src/gflags-download-stamp/gflags-download-patch"
  "deps/gflags/src/gflags-download-stamp/gflags-download-configure"
  "deps/gflags/src/gflags-download-stamp/gflags-download-build"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/gflags-download.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
