########################################################
# build Garfield++ library and examples
# depend on: ROOT
# @author Klaus Zenker
CMAKE_MINIMUM_REQUIRED(VERSION 3.10 FATAL_ERROR)
########################################################

cmake_policy(SET CMP0048 NEW)
PROJECT( Garfield VERSION 0.2.0 LANGUAGES CXX Fortran)

# Say to user he can change the install directory if i didn't do it before
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(STATUS "Install path: ${CMAKE_INSTALL_PREFIX}")
    message(STATUS "If you want to change this path call cmake -DCMAKE_INSTALL_PREFIX=my_install_path ...")
endif()

SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11 CACHE STRING "C++ standard used for compiling" FORCE)
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "C++ standard required" FORCE)
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Allow C++ extensions" FORCE)

### DEPENDENCIES #######################################
IF(DEFINED ENV{ROOTSYS})
    list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
ELSEIF(DEFINED ROOTSYS)
    list(APPEND CMAKE_PREFIX_PATH ${ROOTSYS})
ELSEIF(DEFINED ROOT_CMAKE_DIR)
    list(APPEND CMAKE_PREFIX_PATH ${ROOT_CMAKE_DIR} )
ELSE()
    FIND_PACKAGE( ROOT 5.0 REQUIRED COMPONENTS Geom Gdml QUIET)
    if(NOT ROOT_FOUND)
        MESSAGE(FATAL_ERROR "Please consider to set ROOTSYS or use -DROOTSYS=..."
                        "If still there is a proplem, point to the directory which includes FindROOT.cmake using -DROOT_CMAKE_DIR=...")
    endif()
ENDIF(DEFINED ENV{ROOTSYS})

if(NOT ROOT_FOUND)
FIND_PACKAGE( ROOT 5.0 REQUIRED COMPONENTS Geom Gdml)
endif()

if(ROOT_cxx17_FOUND AND "${CMAKE_CXX_STANDARD}" STREQUAL "11" OR "${CMAKE_CXX_STANDARD}" STREQUAL "14")
	message(WARNING "DQM4hep will be compiled with CXX17")
	set(CMAKE_CXX_STANDARD 17)
elseif(ROOT_cxx14_FOUND AND "${CMAKE_CXX_STANDARD}" STREQUAL "11")
	message(WARNING "DQM4hep will be compiled with CXX14")
	set(CMAKE_CXX_STANDARD 14)
endif()

# Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

# internal root dictionary building works starting with version 6.02
IF(ROOT_VERSION VERSION_LESS 6.02)
    SET( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake )
    INCLUDE( MacroCheckPackageLibs )
    INCLUDE( MacroRootDict )
    MESSAGE(STATUS "Using Garfields FindRoot.cmake")
ENDIF()

## Add root libs ##
MESSAGE(" Root libraries: ${ROOT_LIBRARIES} ")
IF(DEFINED ROOT_INCLUDE_DIRS)
    # ROOT 6
    INCLUDE_DIRECTORIES ( ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} )
ELSE()
    # ROOT 5
    INCLUDE_DIRECTORIES ( ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )
    # COMPONENT mechanism in the FIND_PACKAGE is not supported 
    SET(ROOT_COMPONENT_LIBRARIES "-lGeom -lGdml")
ENDIF(DEFINED ROOT_INCLUDE_DIRS)

## add garfield  #######################################
AUX_SOURCE_DIRECTORY( Source garfield_sources )

## add magboltz ########################################
SET( magboltz_sources Source/magboltz-11.7.f )
SET_SOURCE_FILES_PROPERTIES( ${magboltz_sources}  PROPERTIES LANGUAGE Fortran )

## add Heed ############################################
SET( heed_dirs 
    Heed
    Heed/heed++/code 
    Heed/wcpplib/clhep_units 
    Heed/wcpplib/geometry 
    Heed/wcpplib/ioniz 
    Heed/wcpplib/math 
    Heed/wcpplib/matrix 
    Heed/wcpplib/matter 
    Heed/wcpplib/particle 
    Heed/wcpplib/random 
    Heed/wcpplib/safetl 
    Heed/wcpplib/stream 
    Heed/wcpplib/util )
FOREACH( heed_dir ${heed_dirs} )
    AUX_SOURCE_DIRECTORY( ${heed_dir} heed_sources )
ENDFOREACH()
SET_SOURCE_FILES_PROPERTIES( ${heed_sources} PROPERTIES LANGUAGE CXX)

## Flags to pass to the compiler #######################
SET(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -ansi -Wno-long-long -Woverloaded-virtual -fno-common -pg")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")
SET(CMAKE_CXX_FLAGS_DEBUG  "-ggdb -O0")


ADD_DEFINITIONS( "-c" )


# Add switches used in HEED
# -DNOT_INCLUDE_GASLIB_IN_PACSLIB is used in Heed/heed++/code/PhotoAbsCS.c
SET_SOURCE_FILES_PROPERTIES( ${heed_sources} PROPERTIES COMPILE_FLAGS "-DNOT_INCLUDE_GASLIB_IN_PACSLIB -DGARFIELD_HEED_INTERFACE -DUSE_SRANLUX -DEXCLUDE_FUNCTIONS_WITH_HISTDEF -DINS_CRETURN -DFUNNAMESTACK")

## create dictionary  ##################################
IF(ROOT_VERSION VERSION_LESS 6.02)
    LIST( APPEND ROOT_DICT_INCLUDE_DIRS Include Heed)
    SET( ROOT_DICT_CINT_DEFINITIONS "-DHANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )
    PREPARE_ROOT_DICT_HEADERS( Include )
    GEN_ROOT_DICT_SOURCES( GarfieldDict.cxx )
    LIST( APPEND garfield_sources ${ROOT_DICT_OUTPUT_SOURCES} )
ELSE()
    file (GLOB garfield_header "Include/*.hh")
    ROOT_GENERATE_DICTIONARY( GarfieldDict ${garfield_header} LINKDEF Include/LinkDef.h OPTIONS "-I${PROJECT_SOURCE_DIR}/Heed")
    LIST( APPEND garfield_sources GarfieldDict.cxx )
ENDIF()

## print information on considered files and the root dictionary include directories
MESSAGE( STATUS "###########################################################" )
MESSAGE( STATUS "Considered Garfield files :             ${garfield_sources}" )
MESSAGE( STATUS "-----------------------------------------------------------" )
MESSAGE( STATUS "Considered heed files :                 ${heed_sources}" )
MESSAGE( STATUS "-----------------------------------------------------------" )
MESSAGE( STATUS "Considered magboltz files :             ${magboltz_sources}" )
MESSAGE( STATUS "-----------------------------------------------------------" )
MESSAGE( STATUS "Directory used to find root dict files: ${ROOT_DICT_INCLUDE_DIRS}" )
MESSAGE( STATUS "###########################################################" )


## SET output path for libraries 
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/)
## SET include directories 
INCLUDE_DIRECTORIES( Include Heed )
## build the library
ADD_LIBRARY( ${PROJECT_NAME} SHARED
    ${garfield_sources} 
    ${magboltz_sources}
    ${heed_sources})
    
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES} )
SET_TARGET_PROPERTIES( ${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})

## BUILD examples             ##########################
option(WITH_EXAMPLES "Build Garfield++ examples" OFF)
if(WITH_EXAMPLES)
    MESSAGE(STATUS "Garfield++ exmaples will be build. In order disable this option use -DWITH_EXAMPLES=OFF")
    add_subdirectory(Examples)
else()
    MESSAGE(STATUS "Building of examples is switched off.")
endif()

## EXPORT PACKAGE INFORMATION ##########################
export(PACKAGE ${PROJECT_NAME})

# Configure the config files 
configure_file(CMake/GarfieldConfig.cmake.in
  "${PROJECT_BINARY_DIR}/GarfieldConfig.cmake" @ONLY)
configure_file(CMake/GarfieldConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/GarfieldConfigVersion.cmake" @ONLY)


# decide were the library will end up and to write info in  GarfieldLibraryDepends.cmake
INSTALL(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/ EXPORT GarfieldLibraryDepends)

# decide were to install GarfieldLibraryDepends.camke
INSTALL(EXPORT GarfieldLibraryDepends DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake  COMPONENT cmakefiles)

# now install cofigured config files 
INSTALL( FILES ${PROJECT_BINARY_DIR}/GarfieldConfig.cmake
        ${PROJECT_BINARY_DIR}/GarfieldConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake
        COMPONENT cmakefiles )
        
# install pcm file
file(GLOB PCM_FILE "${PROJECT_BINARY_DIR}/*.pcm")
file(COPY ${PCM_FILE} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

INSTALL (DIRECTORY ${CMAKE_SOURCE_DIR}/Include/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include/Garfield FILES_MATCHING PATTERN "*.h*")
INSTALL (DIRECTORY ${CMAKE_SOURCE_DIR}/Heed/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include/Heed FILES_MATCHING PATTERN "*.h*")
INSTALL (DIRECTORY ${CMAKE_SOURCE_DIR}/Heed/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include/Heed FILES_MATCHING PATTERN "*.ic")
INSTALL (DIRECTORY ${CMAKE_SOURCE_DIR}/Data/ DESTINATION ${CMAKE_INSTALL_PREFIX}/Data)
INSTALL (DIRECTORY ${CMAKE_SOURCE_DIR}/Heed/heed++/database/ DESTINATION ${CMAKE_INSTALL_PREFIX}/Data/Heed)
