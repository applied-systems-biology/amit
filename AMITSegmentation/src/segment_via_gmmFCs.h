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
#include <vector>
#include <opencv2/opencv.hpp>
#include "Region.h"


namespace gmmFCs {

    void fillBR(const std::string outputdir, std::vector<cv::Mat> &images, std::vector<std::vector<Region>> &regions, int &count, const bool &clear_border);
    void classify_PMN(const int &n_classes, std::vector<std::vector<Region>> &binaryRegions, cv::Mat &labels, const std::string &OUTPUTDIR, double &d);
    void create_gmm_image(std::vector<std::vector<Region>> &regions, const int &rows, const int &cols, const cv::Mat &labels, const std::string &outputdir, const int &K, const int &N_TEMPVAR);
    void enhncedOpen (cv::Mat &m, const cv::Mat &kernel);
    void enhncedClose (cv::Mat &m, const cv::Mat &kernel);
    void removeNoiseRegions(cv::Ptr<cv::ml::EM> em_model, std::vector <Region> &regions, const int &rows, const int &cols);
    void combineRegionsOverlap(std::vector<Region> &output, std::vector<Region> &allRegions, std::vector<Region> &classified, std::vector<Region> &morph, std::vector<Region> &added, const double &d);
    void drawImage(cv::Mat & im, std::vector<Region> & regions);
    cv::Mat addingFlat(cv::Ptr<cv::ml::EM> em_model, cv::Mat c, std::vector<Region> binary, std::vector<Region> fungalbinary, std::vector<Region> deadcellbinary, cv::Mat &add, const double &d, const bool &clear_border);
    void perform_gmm_segmentation(const std::vector<cv::Mat> &images, const int &N_THREADS, const int &N_TEMPVAR, const int &C, const double &d, std::vector<std::vector<Region>> binaryRegions, std::vector<std::vector<Region>> fungalbinaryRegions, std::vector<std::vector<Region>> deadcellbinaryRegions, const std::string OUTPUTDIR, const bool &clear_border, const bool &FLAG_DEBUG);

    void gmmFCs_segmentation(const int &n_classes, const std::vector<cv::Mat> &images, std::vector<cv::Mat> &images_fungi, const int &N_THREADS, const int &N_TEMPVAR, const std::string &INPUTDIR, const std::string &OUTPUTDIR, const bool &clear_border, const bool &FLAG_DEBUG);

} // namespace gmmFCs
