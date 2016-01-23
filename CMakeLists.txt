#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ValidateFSData)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
add_definitions(-fpermissive)

# Parallel processing using TOP-C
option(TOPC_USE "Use TOP-C parallel architecture" OFF)
option(TOPC_USE_SEQ "Use TOP-C in sequential mode" OFF)
if(TOPC_USE)
  if(TOPC_USE_SEQ)
    set(PARMODE "--seq")
  else()
    set(PARMODE "--mpi")
  endif()

  add_definitions(-DTOPC_USE=1)
  # Get TOP-C libraries
  execute_process(COMMAND "topc-config" ${PARMODE} "--libs"
          OUTPUT_VARIABLE TOPC_LIBS)

  # Remove trailing whitespace
  string(STRIP ${TOPC_LIBS} TOPC_LIBS)

  # Get TOP-C compiler flags
  execute_process(COMMAND "topc-config" ${PARMODE} "--cflags"
          OUTPUT_VARIABLE TOPC_CFLAGS)
endif()


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc
           ${PROJECT_SOURCE_DIR}/src/*.icc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh
           ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(valFSData ValidateFSData.cc ${sources} ${headers})
target_link_libraries(valFSData ${Geant4_LIBRARIES} ${TOPC_LIBS})


#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(G4_ValidateFSData DEPENDS valFSData)


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS valFSData DESTINATION bin)
