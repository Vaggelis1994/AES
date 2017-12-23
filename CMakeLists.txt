cmake_minimum_required(VERSION 3.6)
project(AES)

set(PROJECT_LINK_LIBS libfplll.so libgmp.so libmpfr.so libqd.so)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lgmp -lfplll -lqd -lmpfr -Wall -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -03 ")
set(SOURCE_FILES main.cpp)

add_executable(AES ${SOURCE_FILES})


target_link_libraries(AES ${PROJECT_LINK_LIBS})
