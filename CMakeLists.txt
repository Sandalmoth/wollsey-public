cmake_minimum_required(VERSION 3.0.0)

set(name wollsey)

add_executable(wollsey main.cpp branching_process.cpp gene.cpp record.cpp protocol.cpp medic.cpp drug.cpp phenotypist.cpp limited_branching_process.cpp constant_branching_process.cpp)
# add_executable(translation test/translation.cpp gene.cpp)
# add_executable(branching test/branching.cpp branching_process.cpp gene.cpp protocol.cpp record.cpp medic.cpp drug.cpp phenotypist.cpp)

# set_target_properties(translation PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test")
# set_target_properties(branching PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test")

# Libraries
# set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5)
if (HDF5_FOUND)
  # TODO don't link when I don't need it
  include_directories(${HDF5_INCLUDE_DIRS})
  target_link_libraries(wollsey ${HDF5_LIBRARIES} -lhdf5_cpp)
  # target_link_libraries(translation ${HDF5_LIBRARIES} -lhdf5_cpp)
  # target_link_libraries(branching ${HDF5_LIBRARIES} -lhdf5_cpp)
endif (HDF5_FOUND)

# Compilation flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++14")

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og -ggdb") # debug compilation
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -Werror")
endif()
