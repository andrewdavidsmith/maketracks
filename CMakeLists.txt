# MIT License
#
# Copyright (c) 2025 Andrew Smith
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.

# to find the version of cmake do
# $ cmake --version
cmake_minimum_required(VERSION 3.30)
project(
  maketracks
  VERSION 0.1.0
  DESCRIPTION
  "Make composite bigWig tracks from transferase methylomes"
  HOMEPAGE_URL https://github.com/andrewdavidsmith/maketracks
  LANGUAGES CXX)

# Set language version used
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_EXTENSIONS off)  # prevents std=gnu++23
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

include(CheckIncludeFileCXX)
include(CheckFunctionExists)
include(CheckCXXCompilerFlag)

include(GNUInstallDirs)
find_package(ZLIB REQUIRED)

# ADS: set the most stringent warnings we can
set(COMPILE_OPTIONS -Wall -Wextra -Wpedantic -Werror -Wfatal-errors)

# Create the CMake target for transferase
add_library(XFR INTERFACE)
set_target_properties(XFR PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES
  "${PROJECT_SOURCE_DIR}/transferase/lib;${PROJECT_SOURCE_DIR}/transferase/cli"
  INTERFACE_LINK_LIBRARIES
  "${PROJECT_SOURCE_DIR}/transferase/build/lib/libtransferase_system.a"
)

add_executable(maketracks maketracks.cpp)
target_compile_options(maketracks PRIVATE ${COMPILE_OPTIONS})
target_link_libraries(maketracks
  PRIVATE
  XFR
  ZLIB::ZLIB
)
