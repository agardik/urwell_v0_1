cmake_minimum_required(VERSION 3.9)
project(urwell_sim)

find_package(Garfield REQUIRED)

add_executable(urwell_sim
  src/main.cc
  src/ProgressBar.cc
  src/DetectorSetup.cc
  src/ElectronSource.cc
  src/AvalancheRunner.cc
  src/Plotting.cc
)

target_include_directories(urwell_sim PRIVATE include)
target_link_libraries(urwell_sim Garfield::Garfield)
