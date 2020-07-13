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


namespace otsu {

    cv::Mat process_otsu(const cv::Mat &src, const int &MIN_REGION_SIZE, const bool &FLAG_DEBUG);
    
    void otsu_segmentation(const std::vector<cv::Mat> &images, const int &N_THREADS, const int &MIN_REGION_SIZE, const std::string &OUT_BINARY, const bool &FLAG_DEBUG);
      
} // namespace otsu
 