#
# Specify list of tests
#

set(testSources
  testConformingVirtualElementOrder1.cpp
   )

set( dependencyList ${parallelDeps} gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  set (dependencyList ${dependencyList} geosx_core)
else()
  set (dependencyList ${dependencyList} ${geosx_core_libs} )
endif()

if (TARGET pugixml::pugixml)
  set( dependencyList ${dependencyList} pugixml::pugixml )
endif()

if (TARGET pugixml)
  set( dependencyList ${dependencyList} pugixml )
endif()

if (TARGET fmt::fmt)
  set( dependencyList ${dependencyList} fmt::fmt )
endif()

if (TARGET fmt)
  set( dependencyList ${dependencyList} fmt )
endif()

#
# Add gtest C++ based tests
#
foreach(test ${testSources})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList}
                      )

    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name}
                  )

    if( ( NOT ${ENABLE_CUDA} ) AND ( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" ) )
        set_target_properties(${test_name} PROPERTIES COMPILE_FLAGS "-Wno-everything")
    endif()

endforeach()

if ( ENABLE_CUDA )
#blt_add_target_compile_flags( TO testShapeFunctionKernel
#                              FLAGS "--maxrregcount=64")
endif()
