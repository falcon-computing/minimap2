FILE(REMOVE_RECURSE
  "CMakeFiles/glog-download"
  "CMakeFiles/glog-download-complete"
  "deps/glog/src/glog-download-stamp/glog-download-install"
  "deps/glog/src/glog-download-stamp/glog-download-mkdir"
  "deps/glog/src/glog-download-stamp/glog-download-download"
  "deps/glog/src/glog-download-stamp/glog-download-update"
  "deps/glog/src/glog-download-stamp/glog-download-patch"
  "deps/glog/src/glog-download-stamp/glog-download-configure"
  "deps/glog/src/glog-download-stamp/glog-download-build"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/glog-download.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
