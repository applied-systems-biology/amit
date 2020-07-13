/*  
Copyright by Jan-Philipp_Praetorius

Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo
Figge
https://www.leibniz-hki.de/en/applied-systems-biology.html
HKI-Center for Systems Biology of Infection
Leibniz Institute for Natural Product Research and Infection Biology -
Hans Knöll Insitute (HKI)
Adolf-Reichwein-Straße 23, 07745 Jena, Germany 

This code is licensed under BSD 2-Clause
See the LICENSE file provided with this code for the full license.
 */

#pragma once

#include <string>
#include <vector>
#include <opencv2/opencv.hpp>


namespace canny {
    
    void gamma_correction(const cv::Mat &src, cv::Mat &dst, const double gamma);
    cv::Mat compute_canny_image(const cv::Mat &src, const bool &FLAG_MINFILTER, const bool &FLAG_MINFILTERPLUSMEDIAN, const int &MORPH_OPEN_CANNY, const int &MORPH_CLOSE_CANNY, const int &MIN_REGION_SIZE, const bool &FLAG_DEBUG);
    
    void canny_segmentation(const std::vector<cv::Mat> &images, const int &N_THREADS, const bool &FLAG_MINFILTER, const bool &FLAG_MINFILTERPLUSMEDIAN, const int &MORPH_OPEN_CANNY, const int &MORPH_CLOSE_CANNY, const int &MIN_REGION_SIZE, const std::string &OUT_GREEN, const bool &FLAG_DEBUG);
    
} // namespace canny
