cmake_minimum_required(VERSION 2.6)
project(discord)

add_executable(discord quit_unitary.cpp su2.cpp opt-func.cpp qubit.cpp quit_jarlskog.cpp quit_import.cpp quit_toolbox.cpp qr_complex.cpp main.cpp)

link_directories(${CMAKE_CURRENT_BINARY_DIR}/lib)


install(TARGETS discord RUNTIME DESTINATION bin)
