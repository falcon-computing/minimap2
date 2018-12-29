message(STATUS "downloading...
     src='https://s3.amazonaws.com/fcs-build-public/falcon-lic-latest.tgz'
     dst='/curr/tianj/minimap2-ori/minimap2-mdbs/build/falconlm/src/falcon-lic-latest.tgz'
     timeout='none'")




file(DOWNLOAD
  "https://s3.amazonaws.com/fcs-build-public/falcon-lic-latest.tgz"
  "/curr/tianj/minimap2-ori/minimap2-mdbs/build/falconlm/src/falcon-lic-latest.tgz"
  SHOW_PROGRESS
  # no EXPECTED_HASH
  # no TIMEOUT
  STATUS status
  LOG log)

list(GET status 0 status_code)
list(GET status 1 status_string)

if(NOT status_code EQUAL 0)
  message(FATAL_ERROR "error: downloading 'https://s3.amazonaws.com/fcs-build-public/falcon-lic-latest.tgz' failed
  status_code: ${status_code}
  status_string: ${status_string}
  log: ${log}
")
endif()

message(STATUS "downloading... done")
