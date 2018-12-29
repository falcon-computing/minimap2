set(file "/curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-falcon.tar.gz")
message(STATUS "verifying file...
     file='${file}'")
set(expect_value "2b1bb4285ef4c8963d5e0e338f1952b8")
file(MD5 "${file}" actual_value)
if("${actual_value}" STREQUAL "${expect_value}")
  message(STATUS "verifying file... done")
else()
  message(FATAL_ERROR "error: MD5 hash of
  ${file}
does not match expected value
  expected: ${expect_value}
    actual: ${actual_value}
")
endif()
