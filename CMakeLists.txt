cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# set(THREADS_PREFER_PTHREAD_FLAG ON)


project(SA VERSION 1.0)

add_executable(SA_run main.cpp)
# find_package(Threads REQUIRED)
# target_link_libraries(GA_run PRIVATE Threads::Threads)