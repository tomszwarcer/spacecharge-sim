# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(IonisationChamber)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(main main.C)
target_link_libraries(main Garfield::Garfield)

add_executable(multigrid multigrid.C)
target_link_libraries(multigrid Garfield::Garfield)

add_executable(timekeeper timekeeper.C)
target_link_libraries(timekeeper Garfield::Garfield)


#---Copy all files locally to the build directory-------------------------------
foreach(_file mesh.mphtxt Potential.txt dielectrics.dat)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()

