cmake_minimum_required(VERSION 3.15)

project(PitchShifterPlugin)

if (WIN32)
  add_definitions(-DNOMINMAX)
endif (WIN32)

set(JUCE_SOURCE_DIR "/JUCE") # set juce source dir
add_subdirectory(${JUCE_SOURCE_DIR} "${CMAKE_BINARY_DIR}/juce")

add_subdirectory(rubberband)
add_subdirectory(plugin)
