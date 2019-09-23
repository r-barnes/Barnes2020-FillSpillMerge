cmake_minimum_required (VERSION 3.9)

project (dephier)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(NetCDF)

add_subdirectory(submodules/richdem)

include_directories(submodules)

add_executable(dephier.exe
  src/main.cpp
)

target_link_libraries(dephier.exe 
  PRIVATE 
    richdem
    netcdf
)

target_compile_features(dephier.exe 
  PUBLIC
    cxx_auto_type
    cxx_std_17
)