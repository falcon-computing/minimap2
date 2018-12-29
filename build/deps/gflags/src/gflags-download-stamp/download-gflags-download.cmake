message(STATUS "downloading...
     src='https://s3.amazonaws.com/fcs-build-public/gflags.tar.gz'
     dst='/curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/gflags/src/gflags.tar.gz'
     timeout='none'")




file(DOWNLOAD
  "https://s3.amazonaws.com/fcs-build-public/gflags.tar.gz"
  "/curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/gflags/src/gflags.tar.gz"
  SHOW_PROGRESS
  EXPECTED_HASH;MD5=1de8187489fbced5cc86c2ba241440e4
  # no TIMEOUT
  STATUS status
  LOG log)

list(GET status 0 status_code)
list(GET status 1 status_string)

if(NOT status_code EQUAL 0)
  message(FATAL_ERROR "error: downloading 'https://s3.amazonaws.com/fcs-build-public/gflags.tar.gz' failed
  status_code: ${status_code}
  status_string: ${status_string}
  log: ${log}
")
endif()

message(STATUS "downloading... done")
