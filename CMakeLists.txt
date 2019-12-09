cmake_minimum_required(VERSION 2.8)

####
# project settings
####

project(MicroLensDemo)

set(sample_VERSION_MAJOR 1)
set(sample_VERSION_MINOR 0)


####
# headers
####

set(HEADERS
)


####
# sources
####

set(SOURCES
	microlens_demo.cpp
)


####
# dependencies
####

find_package(GLAMER NO_MODULE REQUIRED)

include_directories(${GLAMER_INCLUDE_DIRS})


####
# build
####

add_executable(MicroLensDemo ${HEADERS} ${SOURCES})

target_link_libraries(MicroLensDemo ${GLAMER_LIBRARIES})


