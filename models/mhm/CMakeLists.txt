# Upstream's build unchanged, plus HydroTuring's driver linked to its library.
cmake_minimum_required(VERSION 3.18)
project(ht_mhm LANGUAGES Fortran C)
add_subdirectory(/src /src/build/mhm)
add_executable(ht_mhm_driver ht_mhm_driver.f90)
target_link_libraries(ht_mhm_driver PRIVATE mhm_lib)
