add_library(RubberBand RubberBandSingle.cpp)
target_include_directories(RubberBand PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/../rubberband"
        "${CMAKE_CURRENT_SOURCE_DIR}/../src"
        "${CMAKE_CURRENT_SOURCE_DIR}/../src"
        )