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
project(AMITTracking)

set ( CMAKE_C_FLAGS -Werror )
set ( CMAKE_CXX_STANDARD 17)

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif() 

find_package( OpenMP REQUIRED )
find_package( OpenCV REQUIRED core highgui imgproc imgcodecs ml)
find_package( nlohmann_json REQUIRED )
feature_summary(WHAT ALL)

add_executable(AMITTracking
	src/main_tracking.cpp
	src/tracking_region_association.h 
	src/tracking_region_association.cpp
	src/ibp_cluster_detection.h 
	src/ibp_cluster_detection.cpp
)

target_compile_options( AMITTracking PRIVATE -Wall ${OpenMP_CXX_FLAGS})
target_link_libraries( AMITTracking PRIVATE ${OpenMP_CXX_FLAGS})
target_include_directories( AMITTracking PRIVATE ${OpenCV_INCLUDE_DIRS} )
target_link_libraries( AMITTracking PRIVATE OpenMP::OpenMP_CXX ${OpenCV_LIBS} )
target_link_libraries( AMITTracking PUBLIC AMITutils )
