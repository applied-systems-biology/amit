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


