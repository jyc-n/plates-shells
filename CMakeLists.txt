cmake_minimum_required(VERSION 3.9)
project(plates_shells)

include_directories(/usr/local/include)
set(CMAKE_CXX_STANDARD 11)

add_executable(plates_shells main.cpp pre_processor.h pre_processor.cpp parameters.cpp parameters.h geometry.cpp geometry.h solver.cpp solver.h element.cpp element.h node.cpp node.h )