# CMake generated Testfile for 
# Source directory: /curr/tianj/minimap2-ori/minimap2-mdbs/test
# Build directory: /curr/tianj/minimap2-ori/minimap2-mdbs/build/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(minimap-test "bash" "test-driver.sh" "/curr/tianj/minimap2-ori/minimap2-mdbs/build/test/test_minimap")
SET_TESTS_PROPERTIES(minimap-test PROPERTIES  WORKING_DIRECTORY "/curr/tianj/minimap2-ori/minimap2-mdbs/test")
