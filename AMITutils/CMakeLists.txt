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

project(AMITutils)

include(FeatureSummary)

set(CMAKE_CXX_STANDARD 17)

find_package( OpenCV COMPONENTS REQUIRED core highgui imgproc video)
find_package( Boost 1.55.0 COMPONENTS REQUIRED filesystem program_options )
#find_package(LEMON REQUIRED)
feature_summary(WHAT ALL)
message(STATUS "OpenCV libraries: ${OpenCV_LIBS}")
message(STATUS "LEMON libraries: ${LEMON_LIBRARIES}")

add_library(AMITutils
    IOputs.h IOputs.cpp
	Segmentation.h Segmentation.cpp
    Region.h Region.cpp
	RegionSplitting.h RegionSplitting.cpp
	Classification.h Classification.cpp
	Calcs.h Calcs.cpp
	CellTrack.h CellTrack.cpp
	Tracking.h Tracking.cpp
	Outputs.h Outputs.cpp
	RegionP.h RegionP.cpp
	CellTrackP.h CellTrackP.cpp
	ImageProcessingToolbox.h ImageProcessingToolbox.cpp
    visualize.h visualize.cpp
)

target_link_libraries( AMITutils PUBLIC Boost::filesystem Boost::program_options)
target_include_directories( AMITutils INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries( AMITutils PUBLIC ${OpenCV_LIBS})


