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
project(AMIT)

add_subdirectory(AMITutils)
add_subdirectory(AMITSegmentation)
add_subdirectory(AMITTracking)
add_subdirectory(AMITStateClassification)


set(CPACK_GENERATOR "DEB;TGZ")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "ASB")

set(CPACK_DEBIAN_PACKAGE_DEPENDS "libopencv-core (>=3.2)")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS "ON")
set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS "ON")
set(CPACK_PACKAGE_NAME "AMIT")
set(CPACK_PACKAGE_VENDOR "ASB")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "AMIT - algorithm for migration and interaction tracking")
set(CPACK_PACKAGE_VERSION "1.0.0")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")

include(CPack)

