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
#include "RegionP.h"


namespace state_classification
{
    
    double euklideanDistance(cv::Point &v, cv::Point &w, std::optional<double> factor = std::nullopt);
    double computeAverageDistanceFromGreenPixelsToCentroid(RegionP &region);

    void computeFeatures(std::vector<Region> &regions, std::vector<Region> &pmns, std::vector<RegionP> &fungi, cv::Mat &green, cv::Mat &red, std::string &file, const int &t, int ID, const double &f);
    
    std::string stateClassification(std::vector<cv::Mat> images_binary, std::vector<cv::Mat> images_fungi, std::vector<cv::Mat> images_dead, std::string &OUTPURDIR, const int &N_THREADS, const bool &FLAG_DEBUG);

}