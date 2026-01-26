cmake_minimum_required(VERSION 3.16)
project(BQBN-Reborn LANGUAGES CXX)

# -------------------------------
# Set C++ standard
# -------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -------------------------------
# Output directory
# -------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# -------------------------------
# Include directories
# -------------------------------
include_directories(${PROJECT_SOURCE_DIR}/include)

# -------------------------------
# Collect source files
# -------------------------------
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# -------------------------------
# Create executable
# -------------------------------
add_executable(BQBN-Reborn ${SOURCES})


# -------------------------------
# Create updater (Windows only)
# -------------------------------
if(WIN32)
    add_executable(BQBN_updater
        updater/updater.cpp
    )
endif()
