cmake_minimum_required (VERSION 3.4)
project (decovar CXX)

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
# Targets
#--------------------------------------------------------------------------------------------------

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable (decovar main.cpp)
target_link_libraries (decovar sharg::sharg fmt::fmt-header-only biocpp::core biocpp::io)
target_compile_options(decovar PRIVATE -Wall -Wextra)
