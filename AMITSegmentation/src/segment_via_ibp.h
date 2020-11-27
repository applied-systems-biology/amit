/*  
 * Copyright by Jan-Philipp_Praetorius
 * Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
 * https://www.leibniz-hki.de/en/applied-systems-biology.html
 * HKI-Center for Systems Biology of Infection
 * Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
 * Adolf-Reichwein-Straße 23, 07745 Jena, Germany 
 *
 * This code is licensed under BSD 2-Clause
 * See the LICENSE file provided with this code for the full license.
 */

#pragma once

#include <vector>
#include <string>
#include <opencv2/opencv.hpp>


namespace ibp
{

    void houghTransformation(const cv::Mat &src, cv::Mat &dst, const int &minIntersection, const int &minLineLength, const int &maxLineGap, const bool &use_canny, const bool &draw_full_line);
    void myfindContours(const cv::Mat &src, const cv::Mat &original, cv::Mat &dst, const cv::Scalar color = cv::Scalar(255,0,0) );
    void enhance_gridlines(const cv::Mat &src, cv::Mat &dst, const int &min_area_size, const int &minIntersection, const int &minLineLength, const int &maxLineGap);
    void approx_houghLines(const cv::Mat &src, cv::Mat &dst, const bool &FLAG_DEBUG);
    cv::Mat segSplitFull(const cv::Mat &src, const cv::Mat &sd1, const cv::Mat &sd2, const std::string &rmg, const std::string &FFT_MASK_PATH, const bool &clear_border, const int &min_region_size, double &mythresh, const int&erode_kernelSize, const std::string &OUTPUTDIR_CNT, const int &i, const bool &FLAG_DEBUG);

    void ibp_segmentation(const std::vector<cv::Mat> &images, const int &N_THREADS, const std::string &rmg, const std::string &FFT_MASK_PATH, const bool &clear_border, const std::string &OUT_BINARY, const int &min_region_size,  const int&sd1_kernelSize, const int&sd2_kernelSize, double &mythresh, const int&erode_kernelSize, const bool &FLAG_DEBUG);
    
} // namespace ibp