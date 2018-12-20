file(GLOB_RECURSE MM_TEST_LIST src/*.cpp)

add_executable(test_minimap ${MM_TEST_LIST} ${MM_SRC_LIST})

add_dependencies(test_minimap Google FalconLM HTS)


if ( USE_FPGA )
  target_link_libraries(test_minimap
    kflow 
    minimap-c
    ${Boost_LIBRARIES} 
    ${FalconLM_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${ZLIB_LIBRARIES}
    ${CMAKE_DL_LIBS}
    ${HTS_LIBRARIES}
    ${Blaze_LIBRARIES}
    ${KSight_LIBRARIES}
    ${FPGA_LIBRARIES}
    ${Google_LIBRARIES})
else()
  target_link_libraries(test_minimap
    kflow 
    minimap-c
    ${Boost_LIBRARIES} 
    ${FalconLM_LIBRARIES}
    ${Google_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${ZLIB_LIBRARIES}
    ${CMAKE_DL_LIBS}
    ${HTS_LIBRARIES})
endif()

add_test(NAME    minimap-test
         COMMAND bash test-driver.sh $<TARGET_FILE:test_minimap>
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test)
