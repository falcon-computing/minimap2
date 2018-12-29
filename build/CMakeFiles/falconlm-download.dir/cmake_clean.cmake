FILE(REMOVE_RECURSE
  "CMakeFiles/falconlm-download"
  "CMakeFiles/falconlm-download-complete"
  "falconlm/src/falconlm-download-stamp/falconlm-download-install"
  "falconlm/src/falconlm-download-stamp/falconlm-download-mkdir"
  "falconlm/src/falconlm-download-stamp/falconlm-download-download"
  "falconlm/src/falconlm-download-stamp/falconlm-download-update"
  "falconlm/src/falconlm-download-stamp/falconlm-download-patch"
  "falconlm/src/falconlm-download-stamp/falconlm-download-configure"
  "falconlm/src/falconlm-download-stamp/falconlm-download-build"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/falconlm-download.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
