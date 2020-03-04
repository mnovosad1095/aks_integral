cmake_minimum_required(VERSION 3.15)
project(lab3)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
find_package(Threads REQUIRED)
add_executable(lab3 main.cpp integrals.cpp integrals.h)
target_link_libraries(lab3 Threads::Threads)
