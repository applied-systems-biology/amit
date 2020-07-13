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

#include <string>
#include <vector>
#include <boost/optional.hpp>
#include "Region.h"
#include "CellTrack.h"
#include "CellTrackP.h"


namespace outputs {

    /// Interaction Tracking
    void printDataSetFungi(std::string path, std::vector<CellTrack> &tracks);
    void printDataSetPhagocytosis(std::string path, std::vector<CellTrackP> &tracks);
    void showSegmentedRegions2D(std::string &outputdir, std::string &dir_orig_images , int delay, std::vector<std::vector<Region*>> &regions, bool idoutput);
    void showClassifiedSegmentedRegions2D(std::string &outputdir, const std::string &dir_orig_images , const int &delay, std::vector<std::vector<Region>> &regions, const bool &classoutput, const bool &idoutput);

    /// [4.0] NNA tracking
    void showSegmentedRegions2DAndTrajectory(const std::string &outputdir, const std::string &dir_orig_images, const int &delay, int maxt, std::vector<CellTrack> &tracks, const cv::Scalar &trackcolor, const cv::Scalar &idcolor, const bool &showInteractions, const cv::Size2i imagesize);
    /// [5.0] Interaction tracking  &&  [6.0] Combine tracklets globally
    void showSegmentedRegions2DAndTrajectory(const std::string &outputdir, const std::string &dir_orig_images, const int &delay, std::vector<CellTrack> &tracks_single, std::vector<CellTrack> &tracks_cluster);
    void showSegmentedRegions2DAndTrajectory(std::string outputdir, std::string dir_orig_images , int delay,  std::vector<CellTrackP> &tracks, std::vector<CellTrackP> &tracks2);

    /// [8.0] Outputs
    void showSegmentedRegions2DAndTrajectory(std::string &outputdir, std::string dir_orig_images , const int &delay, const int &maxt, std::vector<CellTrack> &tracks, const cv::Scalar &trackcolor, const cv::Scalar &idcolor, const bool showInteractions, cv::Size2i imagesize, const int fontsize);
    void showSegmentedRegions2DAndTrajectory(std::string &outputdir, std::string &dir_orig_images, const int delay, int maxt, std::vector<CellTrack*> &tracks, const cv::Scalar &trackcolor, const cv::Scalar &idcolor, const bool showInteractions, const cv::Size2i imagesize, const int fontsize);

    /// plotTrajectory - overload - functions
    void plotTrajectory(cv::Mat &image, CellTrack &track, const int &t, const cv::Scalar &trackcolor);
    void plotTrajectory(cv::Mat &image, CellTrack *track, const int &t, const cv::Scalar &trackcolor);

    void plotTrajectory(cv::Mat &image, CellTrack &track, const int &t1, const int &t2, const int &tmax, const cv::Scalar &trackcolor);
    void plotTrajectory(cv::Mat &image, CellTrack *track, const int &t1, const int &t2, const int &tmax, const cv::Scalar &trackcolor);
    void plotTrajectory(cv::Mat &image, CellTrackP *track, const int t1, const int t2, const int tmax, const cv::Scalar trackcolor);

    void plotTrajectory(cv::Mat &image, CellTrack &track, const int &t1, const int &t2, const int &tmax);

    void plotTrajectory(cv::Mat &image, CellTrack &track, const int &t);

    /// plotAllTrajectories
    void plotAllTrajectories(const std::string &outfile, std::vector<CellTrackP> &tracks, int cols, int rows);
    void plotAllTrajectories(const std::string &outfile, std::vector<CellTrack> &tracks, const int &cols, const int &rows);

    void printDataSet(const std::string file, std::vector<CellTrack> &tracks);

    void printDataSetRegions(const std::string path, std::vector<CellTrack> &tracks);

    void printClusterRegions(const std::string path, std::vector<CellTrack> &tracks);

} // outputs