message(STATUS "downloading...
     src='https://s3.amazonaws.com/fcs-build-public/glog-falcon.tar.gz'
     dst='/curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-falcon.tar.gz'
     timeout='none'")




file(DOWNLOAD
  "https://s3.amazonaws.com/fcs-build-public/glog-falcon.tar.gz"
  "/curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-falcon.tar.gz"
  SHOW_PROGRESS
  EXPECTED_HASH;MD5=2b1bb4285ef4c8963d5e0e338f1952b8
  # no TIMEOUT
  STATUS status
  LOG log)

list(GET status 0 status_code)
list(GET status 1 status_string)

if(NOT status_code EQUAL 0)
  message(FATAL_ERROR "error: downloading 'https://s3.amazonaws.com/fcs-build-public/glog-falcon.tar.gz' failed
  status_code: ${status_code}
  status_string: ${status_string}
  log: ${log}
")
endif()

message(STATUS "downloading... done")
