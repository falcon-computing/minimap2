FILE(REMOVE_RECURSE
  "CMakeFiles/hts-download"
  "CMakeFiles/hts-download-complete"
  "deps/hts/src/hts-download-stamp/hts-download-install"
  "deps/hts/src/hts-download-stamp/hts-download-mkdir"
  "deps/hts/src/hts-download-stamp/hts-download-download"
  "deps/hts/src/hts-download-stamp/hts-download-update"
  "deps/hts/src/hts-download-stamp/hts-download-patch"
  "deps/hts/src/hts-download-stamp/hts-download-configure"
  "deps/hts/src/hts-download-stamp/hts-download-build"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/hts-download.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
