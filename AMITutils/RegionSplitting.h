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
#include "RegionP.h"


namespace RegionSplitting
{

    void separateClusterGMM(Region &cluster, std::vector<Region> &subregions);
	void separateClusterGMM2(RegionP &cluster, std::vector<RegionP> &subregions);
    
	void associateCellsBipGraphMatching(const int &K, std::vector<Region> &regions, std::vector<Region> &newregions);
    void associateCellsBipGraphMatching2(int K, std::vector<RegionP> &regions, std::vector<RegionP> &newregions);

	void generateCostFunctionPairwiseEuklideanDistance(cv::Mat &cost, std::vector<Region*> & regions1, std::vector<Region*> & regions2);
    void generateCostFunctionPairwiseEuklideanDistanceAndPhagocytosis(cv::Mat &cost, std::vector<RegionP> &regions1, std::vector<RegionP> &regions2);
	
	void minimumWeightedGraphMatching(cv::Mat &cost, std::vector<int> &match);

} // RegionSplitting

template<typename T> void vecToVecP(std::vector<T> & v, std::vector<T*> &res){
	res.clear();
	typename std::vector<T>::iterator it = v.begin();
	for(; it != v.end(); it++){
		res.push_back(&(*it));
	}
}