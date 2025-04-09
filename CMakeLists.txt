# Minimum requirements
cmake_minimum_required (VERSION 3.14)

# Project name
project(Voxeler VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compile program from source files
set(SOURCE_LIST
	"source/main.cpp"

	"source/voxeler.h"
	"source/voxeler.cpp"
)

add_executable(
  voxeler
  ${SOURCE_LIST}
)
target_link_libraries(
  voxeler
)

# Set correct start up project when opening Visual Studio
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT voxeler)