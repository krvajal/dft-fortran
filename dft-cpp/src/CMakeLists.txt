cmake_minimum_required(VERSION 3.5)
project(dft)
# packages
find_package(gsl REQUIRED)
find_package(armadillo REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(LIBS ${LIBS} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})


include_directories(include)
include_directories(odeint/include)
include_directories(armadillo-7.400.2)
include_directories(${GSL_INCLUDE_DIRS})


set(SOURCE_FILES main.cpp )

add_executable(dft ${SOURCE_FILES})
target_link_libraries(dft ${LIBS} ${ARMADILLO_LIBRARIES})
