# basic cmake header
cmake_minimum_required(VERSION 3.1)
project(freesdr)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 14)

#
include_directories(include)
install(DIRECTORY include DESTINATION include)

add_executable(freesdr_test src/test/Waterfall.cpp)
enable_testing()
add_test(NAME Test COMMAND freesdr_test)
