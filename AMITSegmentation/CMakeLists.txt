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
project(AMITSegmentation)

set ( CMAKE_C_FLAGS -Werror )
set ( CMAKE_CXX_STANDARD 17)

find_package( OpenMP REQUIRED )
find_package( OpenCV REQUIRED core highgui imgproc imgcodecs ml)
find_package( nlohmann_json REQUIRED )
feature_summary(WHAT ALL)

add_executable(AMITSegmentation
	src/main_segmentation.cpp
	src/segment_via_GMM.h 
	src/segment_via_GMM.cpp 
	src/segment_via_canny.h 
	src/segment_via_canny.cpp 
	src/segment_via_ibp.h 
	src/segment_via_ibp.cpp 
	src/segment_via_otsu.h 
	src/segment_via_otsu.cpp 
	src/segment_via_gmmFCs.h 
	src/segment_via_gmmFCs.cpp 
	src/segment_via_singleCells.h 
	src/segment_via_singleCells.cpp 
    src/misaxx_convolve.h 
	src/misaxx_convolve.cpp 
)

message(STATUS ${AMITutils})

target_compile_options( AMITSegmentation PRIVATE -Wall ${OpenMP_CXX_FLAGS})
target_link_libraries( AMITSegmentation PRIVATE ${OpenMP_CXX_FLAGS})
target_include_directories( AMITSegmentation PRIVATE ${OpenCV_INCLUDE_DIRS} )
target_link_libraries( AMITSegmentation PRIVATE OpenMP::OpenMP_CXX ${OpenCV_LIBS} )
target_link_libraries( AMITSegmentation PUBLIC AMITutils )

install(TARGETS AMITSegmentation
		RUNTIME
		DESTINATION bin
		COMPONENT application)


