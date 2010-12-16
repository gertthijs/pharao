cmake_minimum_required (VERSION 2.6) 

# Project
project(pharao)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
find_package(OpenBabel2 REQUIRED)

# Version
set (${PROJECT_NAME}_VERSION_MAJOR 3)
set (${PROJECT_NAME}_VERSION_MINOR 0)
set (${PROJECT_NAME}_VERSION_PATCH 3)
configure_file ("${PROJECT_SOURCE_DIR}/pharaoConfig.h.in" "${PROJECT_SOURCE_DIR}/include/pharaoConfig.h")

# Release type
set (CMAKE_BUILD_TYPE Release)

# Sieve code
set (PharaoDir ${PROJECT_SOURCE_DIR}/src)
set (sources     ${PharaoDir}/main.cpp
                 ${PharaoDir}/addBest.cpp
                 ${PharaoDir}/alignment.cpp
                 ${PharaoDir}/aromFuncCalc.cpp
                 ${PharaoDir}/calcPharm.cpp
                 ${PharaoDir}/chargeFuncCalc.cpp
                 ${PharaoDir}/compScore.cpp
                 ${PharaoDir}/coordinate.cpp
                 ${PharaoDir}/functionMapping.cpp
                 ${PharaoDir}/getExt.cpp
                 ${PharaoDir}/hAccFuncCalc.cpp
                 ${PharaoDir}/hDonFuncCalc.cpp
                 ${PharaoDir}/hybridCalc.cpp
                 ${PharaoDir}/lipoFuncCalc.cpp
                 ${PharaoDir}/logOut.cpp
                 ${PharaoDir}/logPharmacophores.cpp
                 ${PharaoDir}/logScores.cpp
                 ${PharaoDir}/mainErr.cpp
                 ${PharaoDir}/mainWar.cpp
                 ${PharaoDir}/options.cpp
                 ${PharaoDir}/parseCommandLine.cpp
                 ${PharaoDir}/pharMerger.cpp
                 ${PharaoDir}/pharmacophore.cpp
                 ${PharaoDir}/printHeader.cpp
                 ${PharaoDir}/printInfo.cpp
                 ${PharaoDir}/printProgress.cpp
                 ${PharaoDir}/printUsage.cpp
                 ${PharaoDir}/result.cpp
                 ${PharaoDir}/siMath.cpp
                 ${PharaoDir}/solutionInfo.cpp
                 ${PharaoDir}/stringTokenizer.cpp
                 ${PharaoDir}/utilities.cpp
)

# Include files from pharao project
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${OPENBABEL2_INCLUDE_DIR})

# Executables
add_executable(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} ${OPENBABEL2_LIBRARIES})
install(TARGETS ${PROJECT_NAME} DESTINATION bin)


#-----------------------------------------
# packaging source tar.gz
#
set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Pharmacophore .")
set(CPACK_PACKAGE_VENDOR "Silicos NV")
set(CPACK_PACKAGE_VERSION_MAJOR @pharao_VERSION_MAJOR@)
set(CPACK_PACKAGE_VERSION_MINOR @pharao_VERSION_MINOR@)
set(CPACK_PACKAGE_VERSION_PATCH @pharao_VERSION_PATCH@)

set(CPACK_RESOURCE_FILE_LICENSE "COPYING")
set(CPACK_RESOURCE_FILE_README "README.TXT")

set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES 
	".svn;/CMakeFiles/;.DS_Store;/_CPack_Packages/;CMakeCache.txt;pharao.tmproj;Makefile"
)
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

include(CPack)

