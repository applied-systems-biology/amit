# /**
# * Copyright by Jan-Philipp Praetorius
# * Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
# * https://www.leibniz-hki.de/en/applied-systems-biology.html
# * HKI-Center for Systems Biology of Infection
# * Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
# * Adolf-Reichwein-Straße 23, 07745 Jena, Germany
# *
# * This code is licensed under BSD 2-Clause
# * See the LICENSE file provided with this code for the full license.
# */

cmake_minimum_required(VERSION 3.10)

set(CMAKE_VERBOSE_MAKEFILE OFF)

include(FeatureSummary)
project(AMITStateClassification)

set ( CMAKE_C_FLAGS -Werror )
set ( CMAKE_CXX_STANDARD 17 )

find_package( OpenMP REQUIRED )
find_package( OpenCV REQUIRED core highgui imgproc imgcodecs ml)
find_package( nlohmann_json REQUIRED )
feature_summary(WHAT ALL)

add_executable(AMITStateClassification
    src/main_state_classification.cpp
    src/state_classification.h
    src/state_classification.cpp
)

message(STATUS ${AMITutils})

target_compile_options( AMITStateClassification PRIVATE -Wall ${OpenMP_CXX_FLAGS})
target_link_libraries( AMITStateClassification PRIVATE ${OpenMP_CXX_FLAGS})
target_include_directories( AMITStateClassification PRIVATE ${OpenCV_INCLUDE_DIRS} )
target_link_libraries( AMITStateClassification PRIVATE OpenMP::OpenMP_CXX ${OpenCV_LIBS} )
target_link_libraries( AMITStateClassification PUBLIC AMITutils )

configure_file(src/classification.R ${CMAKE_CURRENT_BINARY_DIR}/src/classification.R COPYONLY)
configure_file(src/ref7.csv ${CMAKE_CURRENT_BINARY_DIR}/src/ref7.csv COPYONLY)
configure_file(src/data130819_red.csv ${CMAKE_CURRENT_BINARY_DIR}/src/data130819_red.csv COPYONLY)
