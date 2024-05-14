cmake_minimum_required(VERSION 3.27)
project(DA2324_PRJ2_G14_7)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(DA2324_PRJ2_G14_7
        data_structures
        DataReader/DataReader.cpp
        DataReader/DataReader.h
        menu.h
        main.cpp
        menu.cpp
        Class/Edges.cpp)
