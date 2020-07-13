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

#include <vector>
#include <string>
#include "Region.h"


namespace tracking_region_association {

    void setFlatStatus(std::vector<Region> &regions, const bool stat);
    void setEllipsod(std::vector<Region> &regions);

    // estimate Gaussian Mixture distribution of cell areas
    void assign_gmm_class_to_region(std::vector<std::vector<Region>> &regions, const cv::Mat &labels,  const bool &FLAG_DEBUG);
    void estimate_gmm_cell_areas(const int &C, std::vector<std::vector<Region>> &regions, const int &DELAY, const std::string &INPUTDIR_COLOR, const std::string &OUTDIR, const bool &FLAG_DEBUG);

    // assign only single cell class because in the data are no clusters
    void assign_single_class_to_region(std::vector<std::vector<Region>> &regions, const bool &FLAG_DEBUG);

}