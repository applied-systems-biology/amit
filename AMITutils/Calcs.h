/*  
Copyright by Jan-Philipp_Praetorius

Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo
Figge
https://www.leibniz-hki.de/en/applied-systems-biology.html
HKI-Center for Systems Biology of Infection
Leibniz Institute for Natural Product Research and Infection Biology -
Hans Knöll Insitute (HKI)
Adolf-Reichwein-Straße 23, 07745 Jena, Germany 
 */

#pragma once

#include <vector>
#include <opencv2/opencv.hpp>


namespace Calcs {

    void minMaxVec(std::vector<cv::Vec2i> &v, std::vector<cv::Vec2i> &w, cv::Vec2i &min, cv::Vec2i &max);
    
    double computeMean(const std::vector<double> &v);
    double computeVariance(const std::vector<double> &v);
    double computeSd(const std::vector<double> &v);

    bool exists(cv::Vec2i &v, std::vector<cv::Vec2i> &list);

} // Calcs