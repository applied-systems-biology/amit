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


class RegionP : public Region{

	private:
		int I_state; /**< 0 -> single PMN, 1 -> fungus, 2 -> single PMN touching to fungus*/
		bool FC;
		std::vector<int> fungiIDs;
		std::vector<int> deadFungiIDs;
		std::vector<int> phagocytosisIDs;

		cell_class::pathogen klassP;

	public:
		RegionP();
		RegionP(Region &r);
		RegionP(Region &r, bool fc);
		virtual ~RegionP();

		std::vector<cv::Vec2i> green_region_pixels;

		void addFungiId(int id);
		void addDeadFungiId(int id);
		void addPhagocytosisID(int id);
		void addPhagocytosisIDPos(int id, int pos);

		void addGreenPx(const cv::Vec2i &p);
		void addGreenPxs(std::vector<cv::Vec2i> & px);

		void changeFID(int id, int id_new);
		void changePhagocytosisID(int id, int id_new);

		void computeCentroidFromGrayPixels();
		void computeCentroidFromGreenPixels();

		cell_class::pathogen getClassP();
		int getAreaGray();
		cv::Scalar getCentroidScalar();
		cv::Point getCentroidPoint();
		void getDeadFungiIds(std::vector<int> & ids);
		int getIState();
		bool getFCState();
		double getFungalRatioGreen();
		void getFungiIds(std::vector<int> & ids);
		int getNumberOfPhagocytoses();
		void getPhagocytosisIds(std::vector<int> & ids);

		bool hasFungiID(int id);
		bool hasPhagocytozedFungus(int id);

		void removeDeadFungiIds();
		void removeFungiId(int id);
		void removeFungiIds();
		void removeLastPhID();
		void removePhIDs();

		void setClassP(cell_class::pathogen c);
		void setIState(int state);
		void setKlass(cell_class::immune k);
		void setFCState(bool fc);

};
