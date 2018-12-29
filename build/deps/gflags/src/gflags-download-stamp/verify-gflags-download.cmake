set(file "/curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/gflags/src/gflags.tar.gz")
message(STATUS "verifying file...
     file='${file}'")
set(expect_value "1de8187489fbced5cc86c2ba241440e4")
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
