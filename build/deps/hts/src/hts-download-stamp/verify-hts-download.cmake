set(file "/curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/hts/src/htslib-1.3.1.tar.gz")
message(STATUS "verifying file...
     file='${file}'")
set(expect_value "7d0273d4a447b8c01d1a82292eea1ad0")
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
