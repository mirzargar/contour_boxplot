
cmake_minimum_required(VERSION 2.8)

PROJECT(contourboxplot)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(contourboxplot main.cpp Utility.h Utility.cpp CBD.h CBD.cpp CCLabeling.h CCLabeling.cpp ContourExtraction.h ContourExtraction.cpp EnsembleHandler.h EnsembleHandler.cpp Matrix.h SVG.h SVG.cpp)

target_link_libraries(contourboxplot ${ITK_LIBRARIES})
