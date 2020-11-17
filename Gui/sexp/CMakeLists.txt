# SExp - A S-Expression Parser for C++
# Copyright (C) 2015 Ingo Ruhnke <grumbel@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.0)
project(sexp VERSION 0.1.0)

include(GNUInstallDirs)

option(WARNINGS "Switch on extra warnings" OFF)
option(WERROR "Turn warnings into errors" OFF)
option(USE_LOCALE "Use C++ locale support" OFF)
option(USE_CXX17 "Use C++17 features" OFF)

if(WARNINGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(WARNINGS_CXX_FLAGS "${WARNINGS_CXX_FLAGS}"
      -Weverything

      # flags that we deliberately ignore
      -Wno-c++98-compat
      -Wno-c++98-compat-pedantic
      -Wno-float-equal
      -Wno-padded
      -Wno-weak-vtables
      -Wno-disabled-macro-expansion
      -Wno-documentation
      -Wno-reserved-id-macro
      -Wno-sign-conversion
      -Wno-unreachable-code-break
      -Wno-unused-parameter
      -Wno-covered-switch-default
      -Wno-global-constructors
      -Wno-zero-as-null-pointer-constant
      )
  elseif(CMAKE_COMPILER_IS_GNUCXX)
    set(WARNINGS_CXX_FLAGS ${WARNINGS_CXX_FLAGS}
      -pedantic -Wall -Wextra -Wno-c++0x-compat -Wnon-virtual-dtor
      -Wshadow -Wcast-qual -Winit-self -Wno-unused-parameter
      -Weffc++ -Wconversion -Wold-style-cast)
    # -Winline -Wfloat-equal -Wunreachable-code
  endif()
endif()

if(WERROR)
  if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(WARNINGS_CXX_FLAGS ${WARNINGS_CXX_FLAGS}
      -Werror)
  endif()
endif()

if(CMAKE_VERSION VERSION_GREATER 3.6)
  if (USE_CLANG_TIDY)
    find_program(
      CLANG_TIDY_EXE
      NAMES "clang-tidy"
      DOC "Path to clang-tidy executable"
      )

    if(NOT CLANG_TIDY_EXE)
      message(FATAL_ERROR "clang-tidy not found.")
    else()
      message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
      set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
    endif()
  endif()
endif()

add_compile_options(-std=c++11)

if(USE_LOCALE)
  add_definitions(-DSEXP_USE_LOCALE)
endif()

file(GLOB SEXP_SOURCES src/*.cpp)
file(GLOB SEXP_HEADER_SOURCES include/sexp/*.hpp)
add_library(sexp ${SEXP_SOURCES})
set_target_properties(sexp PROPERTIES PUBLIC_HEADER "${SEXP_HEADER_SOURCES}")
target_compile_options(sexp PRIVATE ${WARNINGS_CXX_FLAGS})
target_include_directories(sexp SYSTEM PUBLIC include)

install(TARGETS sexp
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sexp)

if(BUILD_TESTS)
  find_package(Threads REQUIRED)

  # build gtest
  # ${CMAKE_CURRENT_SOURCE_DIR} in include_directories is needed to generate -isystem instead of -I flags
  add_library(gtest_main STATIC ${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googletest/src/gtest_main.cc)
  target_include_directories(gtest_main SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googletest/include/)
  add_library(gtest STATIC ${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googletest/src/gtest-all.cc)
  target_include_directories(gtest SYSTEM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googletest/)
  target_include_directories(gtest SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googletest/include/)
  set_target_properties(gtest_main PROPERTIES CXX_CLANG_TIDY "")
  set_target_properties(gtest PROPERTIES CXX_CLANG_TIDY "")

  # build sexp tests
  file(GLOB TEST_SEXP_SOURCES tests/*.cpp)
  add_executable(test_sexp ${TEST_SEXP_SOURCES})
  target_compile_options(test_sexp PRIVATE ${WARNINGS_CXX_FLAGS})
  target_include_directories(test_sexp PUBLIC src/)
  target_link_libraries(test_sexp
    gtest gtest_main
    sexp
    ${CMAKE_THREAD_LIBS_INIT})

  # add 'make test' target, use 'make test ARGS="-V"' or 'ctest -V' for verbose
  enable_testing()
  add_test(NAME test_sexp
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND test_sexp)
endif()

if(BUILD_BENCHMARKS)
  # build google-benchmark
  include_directories(${PROJECT_SOURCE_DIR}/src)
  file(GLOB BENCHMARK_SOURCES external/benchmark/src/*.cc)
  add_library(benchmark ${BENCHMARK_SOURCES})
  target_compile_definitions(benchmark PRIVATE -DHAVE_STD_REGEX)
  target_compile_options(benchmark PRIVATE -std=c++1y)
  target_include_directories(benchmark SYSTEM PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/external/benchmark/include)
  set_target_properties(benchmark PROPERTIES CXX_CLANG_TIDY "")

  # build benchmarks
  file(GLOB BENCHMARKSOURCES benchmarks/*.cpp)
  foreach(SOURCE ${BENCHMARKSOURCES})
    get_filename_component(SOURCE_BASENAME ${SOURCE} NAME_WE)
    add_executable(${SOURCE_BASENAME} ${SOURCE})
    target_link_libraries(${SOURCE_BASENAME} benchmark sexp ${CMAKE_THREAD_LIBS_INIT})
    set_target_properties(${SOURCE_BASENAME} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/benchmarks/")
    target_compile_options(${SOURCE_BASENAME} PRIVATE -std=c++1y ${WARNINGS_CXX_FLAGS})
  endforeach(SOURCE)
endif()

# EOF #
