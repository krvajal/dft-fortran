cmake_minimum_required(VERSION 3.5)
project(dft)

set(PROJECT_SOURCE_DIR  src)

# packages


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")



include_directories(odeint/include)
include_directories(armadillo-7.400.2/include)



set(SOURCE_FILES src/main.cpp src/verlet.h src/quad/quad.h src/quad/quad_simpson.h)

add_executable(dft ${SOURCE_FILES})

target_link_libraries(dft ${LIBS})
