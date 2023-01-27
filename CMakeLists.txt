# ===========================================================================
#                  deCoVar
# ===========================================================================

cmake_minimum_required (VERSION 3.12.0)
string(ASCII 27 Esc)
set(ColourBold  "${Esc}[1m")
set(ColourReset "${Esc}[m")
set(ColourRed   "${Esc}[31m")

project (decovar CXX)

#--------------------------------------------------------------------------------------------------
# Make "Release" the default cmake build type
#--------------------------------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# Warn if cmake build type is not "Release"
if (NOT CMAKE_BUILD_TYPE STREQUAL Release)
    message (STATUS "${ColourRed}CMAKE_BUILD_TYPE is not \"Release\", your binaries will be slow.${ColourReset}")
endif ()

#--------------------------------------------------------------------------------------------------
# Dependencies
#--------------------------------------------------------------------------------------------------

#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
find_package (biocpp REQUIRED COMPONENTS core io
              HINTS "${CMAKE_SOURCE_DIR}/submodules/biocpp-core/build_system")
find_package (sharg REQUIRED
              HINTS "${CMAKE_SOURCE_DIR}/submodules/sharg-parser/build_system")

set(FMT_DOC OFF)
set(FMT_INSTALL OFF)
set(FMT_TEST OFF)
add_subdirectory(submodules/fmt)

#--------------------------------------------------------------------------------------------------
# Application target
#--------------------------------------------------------------------------------------------------

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable (decovar src/main.cpp src/allele/allele.cpp)
target_link_libraries (decovar sharg::sharg fmt::fmt-header-only biocpp::core biocpp::io)
target_compile_options(decovar PRIVATE -Wall -Wextra)

#--------------------------------------------------------------------------------------------------
# Clang format target
#--------------------------------------------------------------------------------------------------

find_program(CLANG_FORMAT "clang-format-14")

if (CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
    find_program(CLANG_FORMAT "clang-format")
endif()

if (NOT CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
    add_custom_target (check_format "find" "${CMAKE_CURRENT_SOURCE_DIR}/src" "-name" "'*.[ch]pp'" "-exec"
                       ${CLANG_FORMAT} "-style=file" "-n" "-Werror" "{}" "+"
                       COMMENT "Checking the source-code with clang-format.")


    add_custom_target (reformat "find" "${CMAKE_CURRENT_SOURCE_DIR}/src" "-name" "'*.[ch]pp'" "-exec"
                       ${CLANG_FORMAT} "-style=file" "-i" "{}" "+"
                       COMMENT "Reformatting the source-code with clang-format.")
else ()
    message(STATUS "clang_format not found; not adding the 'check_format' and 'reformat' targets.")
endif()
