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
#include "Region.h"


namespace gmm {
    
    cv::Mat create_binary_image_from_regions(const int &rows, const int &cols, std::vector<Region> &regions);
    cv::Mat compute_temporal_variance(const std::vector<cv::Mat> &images);
    cv::Mat compute_local_variance(const cv::Mat &image, const int &kerneldim);
    void classify_pixel(const int &n_i, cv::Mat &res, const std::vector<cv::Mat> &images, cv::Ptr<cv::ml::EM> &em_model, const int &N_TEMPVAR, const int &C, const bool &FLAG_DEBUG);    
    void segment_PMNs_pixel_based_2D(const cv::Mat &image, std::vector<Region> &regions);
    bool is_line(const Region &region);
    void detect_PMNs(const cv::Mat &src, cv::Mat &dst, const int &N_CLOSINGS, const bool &clear_border, const bool &FLAG_DEBUG);
    
    void gmm_segmentation(const std::vector<cv::Mat> &images, const int &N_THREADS, const int &N_TEMPVAR, const int &C, const int &N_CLOSINGS, const std::string &OUT_BINARY, const bool &clear_border, const bool &FLAG_DEBUG);

} // namespace gmm
