cmake_minimum_required(VERSION 3.9)
project(plates_shells)

include_directories(/usr/local/include)
set(CMAKE_CXX_STANDARD 11)

add_executable(plates_shells main.cpp pre_processor.h geometry.cpp geometry.h pre_processor.cpp parameters.cpp parameters.h solver.cpp solver.h stepdata.cpp stepdata.h)