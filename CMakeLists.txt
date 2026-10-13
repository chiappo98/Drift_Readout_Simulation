# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(drift)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

find_package(ROOT REQUIRED COMPONENTS ROOTDataFrame)

# ---Build executable------------------------------------------------------------
add_executable(drift drift.C)
target_link_libraries(drift Garfield::Garfield ROOT::ROOTDataFrame ${ROOT_LIBRARIES})

# add_executable(drift_20x10 drift_20x10.C)
# target_link_libraries(drift_20x10 Garfield::Garfield)

# add_executable(drift_mesh drift_mesh.C)
# target_link_libraries(drift_mesh Garfield::Garfield)

add_executable(generateGas generateGas.C)
target_link_libraries(generateGas Garfield::Garfield)

# add_executable(Fieldlines fieldlines.C)
# target_link_libraries(Fieldlines Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
# foreach(_file ELIST.lis MPLIST.lis NLIST.lis PRNSOL.lis)
#   configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
# endforeach()
