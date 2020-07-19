cmake_minimum_required (VERSION 3.9)

project (dephier)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(NetCDF REQUIRED)

add_subdirectory(submodules/dephier EXCLUDE_FROM_ALL)

add_library(fsm INTERFACE)
target_include_directories(fsm INTERFACE include/)
target_link_libraries(fsm INTERFACE dephier)

include_directories(submodules)

add_executable(fsm_paper_tests.exe paper/fsm_paper_tests.cpp)
target_link_libraries(fsm_paper_tests.exe PRIVATE fsm richdem ${NETCDF_LIBRARIES})
target_compile_features(fsm_paper_tests.exe PUBLIC cxx_auto_type cxx_std_17)

add_executable(fsm_paper_speed_comparison.exe paper/fsm_paper_speed_comparison.cpp)
target_link_libraries(fsm_paper_speed_comparison.exe PRIVATE fsm richdem ${NETCDF_LIBRARIES})
target_compile_features(fsm_paper_speed_comparison.exe PUBLIC cxx_auto_type cxx_std_17)

add_executable(fsm_paper_scaling_tests.exe paper/fsm_paper_scaling_tests.cpp)
target_link_libraries(fsm_paper_scaling_tests.exe PRIVATE fsm richdem ${NETCDF_LIBRARIES})
target_compile_features(fsm_paper_scaling_tests.exe PUBLIC cxx_auto_type cxx_std_17)

add_executable(fsm_unittests.exe
  unittests/test_main.cpp
  unittests/fsm_tests.cpp
)
target_link_libraries(fsm_unittests.exe PRIVATE fsm richdem ${NETCDF_LIBRARIES})
target_compile_features(fsm_unittests.exe PUBLIC cxx_auto_type cxx_std_17)

add_executable(test_temp.exe
  unittests/test_temp.cpp
)
target_link_libraries(test_temp.exe PRIVATE fsm richdem ${NETCDF_LIBRARIES})
target_compile_features(test_temp.exe PUBLIC cxx_auto_type cxx_std_17)


install(TARGETS fsm_paper_tests.exe DESTINATION bin)
install(TARGETS fsm_paper_speed_comparison.exe DESTINATION bin)
