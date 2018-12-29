message(STATUS "downloading...
     src='https://s3.amazonaws.com/fcs-build-public/htslib-1.3.1.tar.gz'
     dst='/curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/hts/src/htslib-1.3.1.tar.gz'
     timeout='none'")




file(DOWNLOAD
  "https://s3.amazonaws.com/fcs-build-public/htslib-1.3.1.tar.gz"
  "/curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/hts/src/htslib-1.3.1.tar.gz"
  SHOW_PROGRESS
  EXPECTED_HASH;MD5=7d0273d4a447b8c01d1a82292eea1ad0
  # no TIMEOUT
  STATUS status
  LOG log)

list(GET status 0 status_code)
list(GET status 1 status_string)

if(NOT status_code EQUAL 0)
  message(FATAL_ERROR "error: downloading 'https://s3.amazonaws.com/fcs-build-public/htslib-1.3.1.tar.gz' failed
  status_code: ${status_code}
  status_string: ${status_string}
  log: ${log}
")
endif()

message(STATUS "downloading... done")
