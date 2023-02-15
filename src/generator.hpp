#pragma once

#include <version>

#if defined(__cpp_lib_generator) && __has_include(<generator>)
#    include <generator>
#else
#    include "../submodules/generator/include/__generator.hpp"
#endif
