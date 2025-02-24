#
# Specify list of tests
#

set( gtest_geosx_tests
     testCompositionalProperties.cpp
     testDamageUtilities.cpp
     testDruckerPrager.cpp
     testElasticIsotropic.cpp
     testKValueInitialization.cpp
     testModifiedCamClay.cpp
     testNegativeTwoPhaseFlash.cpp
     testParticleFluidEnums.cpp
     testPropertyConversions.cpp
     testCubicEOS.cpp
     testRachfordRice.cpp
   )

set( dependencyList gtest constitutive ${parallelDeps} )

if( ENABLE_CUDA_NVTOOLSEXT )
  set( dependencyList ${dependencyList} CUDA::nvToolsExt )
endif()

#
# Add gtest C++ based tests
#
foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
            SOURCES ${test}
            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
            DEPENDS_ON ${dependencyList}
            )

    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name}
                )

endforeach()
