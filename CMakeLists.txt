#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ECRS2021)

#----------------------------------------------------------------------------
# 9/28/2026: Default to an optimized Release build unless the user overrides it.
# A Geant4 application configured with an empty build type is compiled with no
# optimization (-O0), which typically runs several times slower. Build with
# e.g. -DCMAKE_BUILD_TYPE=Debug to override.
#
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Build type: Release, Debug, RelWithDebInfo or MinSizeRel" FORCE)
  message(STATUS "No build type specified - defaulting to Release (optimized).")
endif()


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
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(ECRS2021 ECRS.cc ${sources} ${headers})
# 9/28/2026: The magnetocosmics magnetic-field library path used to be hard-coded
# to a specific user's home directory. Make it an overridable cache variable so
# the project builds on any machine: configure with
#   -DMAGCOS_FIELD_LIB=/path/to/magnetocosmics/fortran/lib/libFortranMagField.a
set(MAGCOS_FIELD_LIB "/home/hexc/G4Work/magnetocosmics/fortran/lib/libFortranMagField.a"
    CACHE FILEPATH "Path to the magnetocosmics libFortranMagField.a")
target_link_libraries(ECRS2021 ${Geant4_LIBRARIES} ${MAGCOS_FIELD_LIB})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B4c. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB4D_SCRIPTS 
  vis.mac
  init_vis.mac
  ECRS_500kRand.txt
  setupECRS_Mag_Field.sh
  )

foreach(_script ${EXAMPLEB4D_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ECRS2021 DESTINATION bin)
