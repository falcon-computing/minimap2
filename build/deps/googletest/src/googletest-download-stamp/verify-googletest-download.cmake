set(file "/curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/googletest/src/googletest.tar.gz")
message(STATUS "verifying file...
     file='${file}'")
set(expect_value "18fda945045354e264e3cca5428525d6")
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
