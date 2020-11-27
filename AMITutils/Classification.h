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

#include <string>
#include <opencv2/opencv.hpp>
#include <vector>
#include "Region.h"


namespace Classification 
{

    void trainSVMonIBPdetection(const std::vector<float> &areas, const std::vector<int> &labels, cv::Ptr<cv::ml::SVM> &model);

    void approxMixDistrEM(const int &n_classes, std::vector<std::vector<Region>> &regions, cv::Mat &labels, const std::string &outdir, bool isFlat);
	void approxMixDistrEM(int K, std::vector<std::vector<Region>> &regions, cv::Mat &labels, cv::Ptr<cv::ml::EM> &em_model, const std::string &outdir, bool isFlat);
    void approxMixDistrEM(int K, std::vector<std::vector<Region*>> &regions, cv::Mat &labels, cv::Ptr<cv::ml::EM> &em_model, const std::string &outdir, bool isFlat);

    void printEMResults(cv::Ptr<cv::ml::EM> &em_model, const std::string &outdir, const std::vector<int> &samplevector, cv::Mat &probs, cv::Mat &likelihood, cv::Mat &labels, bool isFlat);

} // Classification