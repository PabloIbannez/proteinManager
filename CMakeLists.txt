cmake_minimum_required(VERSION 2.4.0)

project(proteinManager)

# Find source and headers files
file(GLOB SOURCES src/*.cpp)
file(GLOB HEADERS include/*.hpp)

# Check libraries
find_package(Boost REQUIRED)

# C++ flags
set(CMAKE_CXX_FLAGS "-fPIC -O3 -std=c++11")

# Include header files
include_directories(include)

# Create shared library
add_library(${PROJECT_NAME} SHARED ${SOURCES})

# Install library
install(TARGETS ${PROJECT_NAME} DESTINATION /usr/lib)

# Install library headers
install(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME})
