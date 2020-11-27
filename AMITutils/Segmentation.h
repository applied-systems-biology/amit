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
#include <opencv2/opencv.hpp>
#include "Region.h"


namespace Segmentation
{

    void getRegion2DIterative(cv::Mat &im, int i, int j, std::vector<cv::Vec2i> &pixels);
    void getRegion2DRecursive(cv::Mat &im, int i, int j, std::vector<cv::Vec2i> &pixels);
    void getContour2D(cv::Mat &im, const std::vector<cv::Vec2i> &region_pixels, std::vector<cv::Vec2i> &pixels);
    bool touchesBorder(const std::vector<cv::Vec2i> &pixels, const cv::Mat &image);
    void fill_holes(std::vector<Region> &regions, const bool &clear_border);

    void segmentRegionPixelBased2DIterative(const cv::Mat &image, std::vector<Region> &regions, const bool &clear_border);
    void segmentRegionPixelBased2DRecursive(const cv::Mat &image, std::vector<Region> &regions);
    
} // Segmentation
