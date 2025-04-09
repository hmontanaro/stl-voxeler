# Minimum requirements
cmake_minimum_required (VERSION 3.14)

# Project name
project(Voxeler VERSION 1.0.0 LANGUAGES C CXX)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE  
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

# Compile program from source files
set(SOURCE_LIST
  "source/main.cpp"

  "source/voxeler.h"
)

add_library(voxeler_lib STATIC
  "source/voxeler.cpp"
)

add_executable(
  voxeler
  ${SOURCE_LIST}
)
target_link_libraries(
  voxeler PRIVATE
  voxeler_lib
)

# Compile program from source files
set(TEST_SOURCE_LIST
  "source/tests/voxeler_test.cpp"
)

add_executable(
  voxeler_test
  ${TEST_SOURCE_LIST}
)
target_link_libraries(
  voxeler_test
  GTest::gtest_main
  voxeler_lib
)

include(GoogleTest)
gtest_discover_tests(voxeler_test)

# Set correct start up project when opening Visual Studio
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT voxeler)