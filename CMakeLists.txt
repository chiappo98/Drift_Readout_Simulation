# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(drift)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

find_package(ROOT REQUIRED COMPONENTS ROOTDataFrame)
find_package(nlohmann_json 3.2.0 CONFIG REQUIRED)

# ---Build executable------------------------------------------------------------
add_executable(drift drift.C)
target_link_libraries(drift Garfield::Garfield ROOT::ROOTDataFrame nlohmann_json::nlohmann_json ${ROOT_LIBRARIES})

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/config")
configure_file(config/config.json "${CMAKE_CURRENT_BINARY_DIR}/config/config.json" COPYONLY)
configure_file(config/simulation.json "${CMAKE_CURRENT_BINARY_DIR}/config/simulation.json" COPYONLY)

add_executable(generateGas generateGas.C)
target_link_libraries(generateGas Garfield::Garfield)

# add_executable(Fieldlines fieldlines.C)
# target_link_libraries(Fieldlines Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
# foreach(_file ELIST.lis MPLIST.lis NLIST.lis PRNSOL.lis)
#   configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
# endforeach()
