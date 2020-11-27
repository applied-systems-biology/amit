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
#include <boost/optional.hpp>
#include "Region.h"
#include "CellTrack.h"
#include "CellTrackP.h"


float const minFloat = -std::numeric_limits<float>::max();

namespace Tracking {

    // [4.0] NNA tracking
    void findNext(Region &current, std::vector<std::vector<Region>> &regions, int t, int id, CellTrack &track);
    
    void trackWithOverlap(std::vector<std::vector<Region>> &regions, std::vector<CellTrack> &tracks);
    void trackWithOverlap(std::vector<std::vector<RegionP>> &regions, std::vector<cv::Mat> images, std::vector<CellTrackP> &tracks);
    void trackWithOverlapClass2(std::vector<std::vector<Region>> &regions, std::vector<CellTrack> &tracks);

    double getAreaDeviation(std::vector<CellTrack> &tracks1);
    double getAreaDeviation(std::vector<CellTrackP> &tracks1);
    double getMaxSpeed(std::vector<CellTrack> &tracks1);
    double getMaxSpeed(std::vector<CellTrackP> &tracks1);
    
    int getTracksLengthTotal(std::vector<CellTrack> &tracks);
    int getTracksLengthTotal(std::vector<CellTrack*> &tracks);

    // [5.0] Interaction tracking for "trackInteractingRegions"
    void addRegionToTrack(Region &r, std::vector<CellTrack> &ct, const int &id, const int &t = -1);
    void addRegionToTrack(RegionP &r, std::vector<CellTrackP> &ct, int id, int t);
    void addInteractionToTrack(std::vector<CellTrack> & tracks, const int trackid, const int t, std::vector<Region> &interacting_cells);
    void addInteractionToTrack(std::vector<CellTrackP> &tracks, const int trackid, const int t, std::vector<RegionP> &interacting_cells);    
    void deleteRegion(std::vector<std::vector<Region>> &regions, const int &t, const int &ID);
    void deleteRegion(std::vector<std::vector<RegionP>> &regions, int t, int ID);
    void reduceOverlaps(Region &cluster, std::vector<Region> &overlaps);
    void reduceOverlaps(Region *cluster, std::vector<Region*> &overlaps);

    void computeOverlapRegionsForCellTracks(std::vector<CellTrack> &tracks, std::vector<CellTrack> &tracks_class2, int t_max);
    void computeOverlapRegionsForCellTracks(std::vector<CellTrackP> &tracks, std::vector<CellTrackP> &tracks_class2, int t_max);
    void computeOverlapRegionsForCellTracksFC(std::vector<CellTrackP> &tracks_fc, std::vector<CellTrackP> &tracks_class2_fc, std::vector<CellTrackP> &tracks_class2, int t_max);
    void computeOverlapRegionsForCellTrack(CellTrackP &track, std::vector<CellTrackP> &tracks, std::vector<CellTrackP> &tracks_class2, int t_max);
	
    bool splitClusterDependentOnArea(int t, double areadev, double c, std::vector<Region> &overlap, std::vector<CellTrack> &tracks, CellTrack &track2, std::vector<std::vector<Region>> &regions);

    void trackInteractingRegions(std::vector<std::vector<Region>> &regions, std::vector<CellTrack> &tracks, std::vector<CellTrack> &tracks_class2, double areadev, double c);
    void createKlassImages(std::vector<cv::Mat> &images_klass0, std::vector<cv::Mat> &images_klass1, std::vector<cv::Mat> &images_klass2, std::vector<std::vector<RegionP>> &regions_over_time, int cols, int rows);
    void drawRegionOnImage(cv::Mat &image, const Region &region);

    // for "combineCellTracks2" 
    void combineCellTracksWithGapSize(std::vector<CellTrack> &tracks1, const int t_max, const int gapsize, const int maxdistance, int rows, int cols, const bool &DEBUG);
    
    void combineTracksWithGraph(std::vector<CellTrack> &tracks1, std::vector<int> &id_start, std::vector<int> &id_end, std::vector<Region> &reg_start, std::vector<Region> &reg_end, std::vector<int> gapsizes, double maxspeed, int rows, int cols, const bool &DEBUG);
    void combineTracksWithGraph(std::vector<CellTrackP> &tracks1, std::vector<int> & id_start, std::vector<int> &id_end, std::vector<RegionP> &reg_start, std::vector<RegionP> &reg_end, std::vector<int> gapsizes, double maxspeed);
    void combineTracksWithGraph2(std::vector<CellTrackP> &tracks1, std::vector<int> & id_start, std::vector<int> & id_end, std::vector<RegionP> & reg_start, std::vector<RegionP> & reg_end, std::vector<int> gapsizes, double maxspeed);
    void generateCostFunctionPairwiseEuklideanDistanceAndPhagocytosis(cv::Mat &cost, std::vector<RegionP> &regions1, std::vector<RegionP> &regions2);

    void combineCellTracks(std::vector<CellTrackP> &tracks1, const int t_max, const int maxgapsize);
    void combineCellTracks2(std::vector<CellTrackP> &tracks1, const int t_max, const int maxgapsize);

    // for "combineCellTracks"    
    void combineTracksWithGraph(std::vector<CellTrack> &tracks1, std::vector<int> &id_start, std::vector<int> &id_end, std::vector<Region> &reg_start, std::vector<Region> &reg_end, std::vector<int> &gapsizes, const double &maxspeed);
    void combineTwoCellTracks(std::vector<CellTrack> &tracks1, const int &id1, const int &id2);
    void combineTwoCellTracks(std::vector<CellTrackP> &tracks1, int id1, int id2);
    
    void computeDistances(cv::Mat & dist, const std::vector<int> &id_start, const std::vector<int> &id_end, std::vector<Region> &reg_start, std::vector<Region> &reg_end);        
    void computeDistances(cv::Mat & dist, const std::vector<int> &id_start, const std::vector<int> &id_end, std::vector<RegionP> &reg_start, std::vector<RegionP> &reg_end);
    void computeDistancesToImageBoundaries(cv::Mat &dist, const std::vector<int> &id_start, const std::vector<int> &id_end, std::vector<Region> &reg_start, std::vector<Region> &reg_end, int rows, int cols);
    double computeDistanceToImageBoundaries(Region &r, Region &p, int rows, int cols);
    
    void combineCellTracks(std::vector<CellTrack> &tracks1,  const int &t_max, const int &maxgapsize);   
    void combineCellTracks(std::vector<CellTrack> &tracks1, const int t_max, const int maxgapsize, double maxspeed, int rows, int cols, const bool &DEBUG);

    void postProcessingClusterSplitting(std::vector<std::vector<Region>> &regions, std::vector<CellTrack> &tracks_singlecells, std::vector<CellTrack> &tracks_cluster);

    // [6.0] Combine tracklets globally
    void combineCellTracks(std::vector<CellTrack> &tracks1, const int &t_min, const int &t_max, const int &maxgapsize);

    // [7.0] Post-processing
    void postProcessing(std::vector<CellTrack> &tracks, const int &maxgapsize);
    void postProcessing(std::vector<CellTrackP> &tracks, int maxgapsize);
    void correctIds(std::vector<CellTrack> &tracks);
    int getNumberOfInteractions(std::vector<CellTrack> &tracks);

} // Tracking