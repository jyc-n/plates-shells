# version requirement
cmake_minimum_required(VERSION 3.9)

# project name
project(plates_shells)

# all source codes
aux_source_directory(${PROJECT_SOURCE_DIR}/src srcs)

# include directories, and header files
include_directories(/usr/local/include)
include_directories(${PROJECT_SOURCE_DIR}/include)

# compiler settings
set(CMAKE_CXX_STANDARD 11)

if(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
    message("Building on Mac OS")
elseif(${CMAKE_SYSTEM_NAME} MATCHES Linux)
    message("Building on Linux Ubuntu")
endif()

# generate executable file
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
add_executable(plates_shells ${srcs})

# link libraries