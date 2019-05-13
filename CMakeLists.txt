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
set(CMAKE_CXX_STANDARD 14)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

# generate executable file
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
add_executable(plates_shells ${srcs})

# Mac OS
if(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
    message("=====Building on Mac OS=====")
    # lapack
    find_package(LAPACK REQUIRED)
    find_library(
        PARADISO_LIB
        NAMES libpardiso600-MACOS-X86-64 pardiso600-MACOS-X86-64
        HINTS "${PROJECT_SOURCE_DIR}/lib"
        NO_DEFAULT_PATH
    )
    message(STATUS "PARADISO_LIB: [${PARADISO_LIB}]")
    find_library(
        OMP_LIB
        NAMES libiomp5 iomp5
        HINTS "${PROJECT_SOURCE_DIR}/lib"
        NO_DEFAULT_PATH
    )
    message(STATUS "OMP_LIB: [${OMP_LIB}]")
    # link libraries
    target_link_libraries(plates_shells ${PARADISO_LIB} ${LAPACK_LIBRARIE} ${OMP_LIB})

# Linux Ubuntu
elseif(${CMAKE_SYSTEM_NAME} MATCHES Linux)
    message("=====Building on Linux Ubuntu=====")
    # linking directories
    link_directories(
                    /usr/lib/x86_64-linux-gnu
                    /usr/lib/x86_64-linux-gnu/lapack
                    /usr/lib/gcc/x86_64-linux-gnu/7.3.0)
    # lapack
    find_package(LAPACK REQUIRED)
    # OpenMP
    find_package(OpenMP REQUIRED)
    if (OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
    # POXIS thread
    find_package(Threads REQUIRED)
    set(pthread ${CMAKE_THREAD_LIBS_INIT})

    target_link_libraries(plates_shells pardiso600-GNU720-X86-64 ${LAPACK_LIBRARIES} gfortran pthread)
endif()