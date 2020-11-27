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

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <nlohmann/json.hpp>
#include <cassert>
#include "IOputs.h"
#include "Region.h"
#include "Segmentation.h"
#include "CellTrack.h"
#include "Tracking.h"
#include "tracking_region_association.h"
#include "ibp_cluster_detection.h"
#include "Outputs.h"

using std::chrono::system_clock;
using std::chrono::duration;
using json = nlohmann::json; 
namespace tra = tracking_region_association;
#define assertm(exp, msg) assert(((void)msg, exp))


auto start = system_clock::now();


/************************************************************/
///                                                          /
///                           main                           /
///                                                          /
/************************************************************/
int main(int argc, char* argv[]) {

    /******************************************************************************/
    //            parse input arguments of command line and get JSON               /
    /******************************************************************************/
    json J = io::get_JSON(argc, argv);

    /// print information, like date, working directory, start time measurement
    time_t rawtime;
  	struct tm *timeinfo;
  	time (&rawtime);
  	timeinfo = localtime(&rawtime);    

    std::cout << " [0.0] Start Tracking at Date:\t\t" << asctime(timeinfo) << std::endl;
    std::cout << " [0.1] Working directory:\t\t" << get_current_dir_name() << std::endl;

    /******************************************************************************/
    ///      read input images, other parameters, create output directory          /
    /******************************************************************************/

    std::string INPUTDIR_ORIGINAL = J["input_original"].get<std::string>();
    std::string INPUTDIR_BINARY = J["input_binary"].get<std::string>();
    std::string OUTPUTDIR = J["output"].get<std::string>();
        
    bool CLEAR_BORDER = J["clear_border"].get<bool>();
    std::string CLUSTER_DETECTION_METHOD = J["cluster_detection_method"].get<std::string>();
    int C = J["number_classes"].get<int>();
    float overlap_thr = J["overlap_threshold"].get<float>();
    int MAXGAPSIZE = J["gap"].get<int>();
    int DELAY = J["delay"].get<int>();    

    size_t N_THREADS = J["n_threads"].get<int>();
    bool FLAG_DEBUG = J["debug"].get<bool>();   
    
    assertm( CLUSTER_DETECTION_METHOD == "ibp" || CLUSTER_DETECTION_METHOD == "gmm" || CLUSTER_DETECTION_METHOD == "none", "Please choose either <ibp> , <gmm> or <none> for cluster detection method");
    if( CLUSTER_DETECTION_METHOD != "ibp" && CLUSTER_DETECTION_METHOD != "gmm" && CLUSTER_DETECTION_METHOD != "none" ){
        std::cout << " [ERROR] Please choose either <ibp> , <gmm> or <none> for parameter cluster_detection_method!" << std::endl;
        return -1;
    }

    std::vector<cv::Mat> images_original, images_binary;

    io::read_images(INPUTDIR_ORIGINAL, images_original, cv::IMREAD_GRAYSCALE);
    io::read_images(INPUTDIR_BINARY, images_binary, cv::IMREAD_GRAYSCALE);

    std::cout << " [0.2.1] Number original-images:\t" << images_original.size() << std::endl;
    std::cout << " [0.2.2] Number binary-images:\t\t" << images_binary.size() << std::endl;

    (void) io::create_directory(OUTPUTDIR);

    /******************************************************************************/
    ///                            segment all regions                             /
    /******************************************************************************/
    std::cout << "\n [2.0] Segmentation of regions" << std::endl; 
	std::vector< std::vector <Region> > regions_over_time;

    /// initialize ROWS and COLS
    int ROWS, COLS;

    if( images_binary.size() > 0 ){
	    ROWS =  images_binary[0].rows;
	 	COLS = images_binary[0].cols;
    }
    else {
        std::cout << " [2.1] ERROR: no binary-images found" << std::endl;
	    return 1;
    }

    /// segmentation
    int n_regions = 0;

	std::vector<std::string>::iterator it, it2;

    for (auto img_gray : images_binary){
        /// get regions via segmentation with regions of current image
		std::vector<Region> regions;
		Segmentation::segmentRegionPixelBased2DIterative(img_gray, regions, CLEAR_BORDER);
        
        tra::setFlatStatus(regions, false);
		tra::setEllipsod(regions);
		
        regions_over_time.push_back(regions);
		n_regions += regions.size();
    }

    std::cout << " [2.1] number of segmented regions:\t" << n_regions << std::endl;

    /******************************************************************************/
    ///  assign cells to single cells or cluster via {gmm, ibp_cluster_detection}  /
    /******************************************************************************/
    std::cout << "\n [3.0] Region classification via:\t" << CLUSTER_DETECTION_METHOD << std::endl;

    /*** estimate Gaussian Mixture distribution of cell areas ***/
    /// classes: (noise) , single cells , cluster 
    if(CLUSTER_DETECTION_METHOD == "gmm"){
        
        tra::estimate_gmm_cell_areas(C, regions_over_time, DELAY, INPUTDIR_ORIGINAL, OUTPUTDIR, FLAG_DEBUG);

    }
    /***   detect clusters and retrieve information from objList.type for each region   ***/
    else if(CLUSTER_DETECTION_METHOD == "ibp"){

        std::vector<cv::Mat> images_objMaps;
        feature_data objList;

        std::cout << " [3.0.0] Overlap-fraction of ibp:\t" << overlap_thr << std::endl;
        ibp_cluster_detection::clustDetectRun(images_binary, images_objMaps, objList, overlap_thr, FLAG_DEBUG);

        /// assign detected cluster / single-cells to the corresponding regions
        ibp_cluster_detection::assign_ibp_class_to_region(regions_over_time, images_objMaps, DELAY, INPUTDIR_ORIGINAL, OUTPUTDIR, FLAG_DEBUG);

        /// write cluster regions to region-files if debug mode is switched ON
        if (FLAG_DEBUG){
            std::string file_regions = io::create_directory(OUTPUTDIR, "region_files/");

            std::vector<CellTrack> tmp_tracks_class_cluster;

            Tracking::trackWithOverlapClass2(regions_over_time, tmp_tracks_class_cluster);

            outputs::printClusterRegions(file_regions, tmp_tracks_class_cluster);
            
            std::cout << " [3.0.1] stored printDataSetRegions" << std::endl;
        }

        /// split the cluster via an iterative watershed segmentation
        std::vector<cv::Mat> images_cluster_splitted;

        int canny_threshold = J["canny_threshold"];

        std::string path_cS = OUTPUTDIR + "/3_ibp_clusterSplitting/";
        ibp_cluster_detection::cluster_splitting(images_original, images_binary, objList, images_cluster_splitted, canny_threshold, path_cS, N_THREADS, FLAG_DEBUG);

        std::cout << " [3.1] cluster splitting done" << std::endl;

        /// redo the regions detection and segmentation
        regions_over_time.clear();

        n_regions = 0;

        for (const auto& img_gray : images_cluster_splitted){
            /// get regions via segmentation with regions of current image
            std::vector<Region> regions;
            Segmentation::segmentRegionPixelBased2DIterative(img_gray, regions, CLEAR_BORDER);

            tra::setFlatStatus(regions, false);
            tra::setEllipsod(regions);

            regions_over_time.push_back(regions);
            n_regions += regions.size();
        }

        std::cout << " [3.2] redo - number of segmented regions:\t" << n_regions << std::endl;

    }
    /*** no clusters are in the underlying binary images (logical not possible or removed during pre-processing) ***/
    else if(CLUSTER_DETECTION_METHOD == "none"){
        
        tra::assign_single_class_to_region(regions_over_time, FLAG_DEBUG);

    }
    else{
        std::cout << " [3.1] ERROR: choose one of the following cluster detection methods: { gmm, ibp, none } " << std::endl;
        return 1;
    }

    /******************************************************************************/
    ///                  track cells that can be tracked for sure                  /
    ///        (only cells that belong to class -> single cells)                   /
    /******************************************************************************/
    std::cout << "\n [4.0] NNA tracking" << std::endl;
	
    std::vector<CellTrack> tracks_class_singlecells;
    std::vector<CellTrack> tracks_class_cluster;
    
    Tracking::trackWithOverlap(regions_over_time, tracks_class_singlecells);
    Tracking::trackWithOverlapClass2(regions_over_time, tracks_class_cluster);

	std::cout << " [4.1] Cell tracks type single cells:\t\t\t" << tracks_class_singlecells.size() << std::endl;
	std::cout << " [4.2] Cell tracks type cell cluster:\t\t\t" << tracks_class_cluster.size() << std::endl;
        
    if(FLAG_DEBUG){
        std::cout << " [4.2.2] output: NNA tracking single cells" << std::endl;
		
        std::string file2 = OUTPUTDIR + "/4_NNA_tracking/";
        cv::Scalar track_color = cv::Scalar(255,255,255);
        cv::Scalar id_color = cv::Scalar(0,0,255);

     	outputs::showSegmentedRegions2DAndTrajectory(file2, INPUTDIR_ORIGINAL, DELAY , 0, tracks_class_singlecells, track_color, id_color, true, cv::Size2i());
	}

    /******************************************************************************/
    ///        manage interaction tracking between single cells and cluster        /
    /******************************************************************************/
    std::cout << "\n [5.0] Interaction tracking" << std::endl;

    /// compute area sd of single cells
	double areadev =  Tracking::getAreaDeviation(tracks_class_singlecells);

	if(FLAG_DEBUG){
		std::cout << " [5.0.1] Std of all cell track areas:\t\t" << areadev << std::endl;
    }

	double mSpeed = Tracking::getMaxSpeed(tracks_class_singlecells);
	std::cout << " [5.1] 1. Speed of cells:\t\t\t"<< mSpeed << std::endl;

    int iterations_w = 0;
    bool change;

    /// check for non - cluster tracks
    if ( Tracking::getTracksLengthTotal(tracks_class_cluster) == 0 )
        change = false;
    else
        change = true;


    while(change){
        /// get total track length of cell clusters
        int N2 = Tracking::getTracksLengthTotal(tracks_class_cluster);

        std::cout << "\n [5.2] Iteration:\t\t\t\t" << iterations_w++ << std::endl;

        /// split the cluster
        Tracking::trackInteractingRegions(regions_over_time, tracks_class_singlecells, tracks_class_cluster, areadev, 2.);

        Tracking::combineCellTracks(tracks_class_singlecells, regions_over_time.size() - 1 , MAXGAPSIZE);

        if(Tracking::getTracksLengthTotal(tracks_class_cluster) == N2){
            change = false;
        }

        std::cout << " [5.2.1] Cell tracks type single cells:\t\t" << tracks_class_singlecells.size() << std::endl;
        std::cout << " [5.2.2] Cell tracks type cell cluster:\t\t" << tracks_class_cluster.size() << std::endl;
    }

    if(FLAG_DEBUG){
        std::string file_nc =  OUTPUTDIR + "/5_tracking_not_combined/";

        outputs::showSegmentedRegions2DAndTrajectory(file_nc, INPUTDIR_ORIGINAL , DELAY, tracks_class_singlecells, tracks_class_cluster);
    }

    std::cout << " [5.3] 2. Max speed is:\t\t\t\t"<< Tracking::getMaxSpeed(tracks_class_singlecells)<< std::endl;
	    
    std::cout << " [5.4] Perform postprocessing" << std::endl;
	
    Tracking::postProcessingClusterSplitting(regions_over_time, tracks_class_singlecells, tracks_class_cluster);

    std::cout << " [5.5] 3. Max speed is:\t\t\t\t"<< Tracking::getMaxSpeed(tracks_class_singlecells)<< std::endl;
	
    /******************************************************************************/
    ///                      Combine tracklets globally                            /
    /******************************************************************************/
    std::cout << "\n [6.0] Combine tracklets globally" << std::endl;

    Tracking::combineCellTracks(tracks_class_singlecells, regions_over_time.size()-1, MAXGAPSIZE, mSpeed);

    std::cout << " [6.1] Cell tracks type single cells:\t\t" << tracks_class_singlecells.size() << std::endl;
	std::cout << " [6.2] Cell tracks type cell cluster:\t\t" << tracks_class_cluster.size() << std::endl;

 	if(FLAG_DEBUG){
 		std::string out = OUTPUTDIR + "/6_tracklet_matching/";

 		outputs::showSegmentedRegions2DAndTrajectory(out, INPUTDIR_ORIGINAL , DELAY, tracks_class_singlecells, tracks_class_cluster);
 	}

    /******************************************************************************/
    ///                                  Postprocessing                            /
    /******************************************************************************/
	std::cout << "\n [7.0] Post-processing " << std::endl;

	Tracking::postProcessing(tracks_class_singlecells, MAXGAPSIZE);
	
	Tracking::correctIds(tracks_class_singlecells);

    std::cout << " [7.1] Cell tracks type single cells:\t\t" << tracks_class_singlecells.size() << std::endl;
	std::cout << " [7.2] Cell tracks type cell cluster:\t\t" << tracks_class_cluster.size() << std::endl;

	int n_interactions = Tracking::getNumberOfInteractions(tracks_class_singlecells);
	
	std::cout << " [7.3] Number of Interactions:\t\t\t" << n_interactions << std::endl;

	/// remove too short tracks
	for(size_t j = 0; j < tracks_class_singlecells.size();j++){

	    if(tracks_class_singlecells.at(j).getLength() <= MAXGAPSIZE ){
	 	    tracks_class_singlecells.erase( tracks_class_singlecells.begin()+j );
	 		j--;
	 	}

	}
    
	std::cout << " [7.4] 4. Max speed is:\t\t\t\t"<< Tracking::getMaxSpeed(tracks_class_singlecells) << std::endl;

    /******************************************************************************/
    ///                                    Outputs                                 /
    /******************************************************************************/
    std::cout << "\n [8.0] Outputs" << std::endl;

    std::string orig_out = OUTPUTDIR + "/tracked/";

    cv::Scalar track_color = cv::Scalar(255,255,255);
    cv::Scalar id_color = cv::Scalar(0,0,255);
    int FONT_SIZE = J["fontsize"].get<int>();
    int MAX_T = 0;
    
	outputs::showSegmentedRegions2DAndTrajectory(orig_out, INPUTDIR_ORIGINAL, DELAY , MAX_T, tracks_class_singlecells, track_color, id_color, true, cv::Size2i(), FONT_SIZE);

	std::cout << " [8.1] stored showSegmentedRegions2DAndTrajectory" << std::endl;

    outputs::plotAllTrajectories(OUTPUTDIR+"/trajectories.png", tracks_class_singlecells, COLS, ROWS);

    std::cout << " [8.2] stored plotAllTrajectories" << std::endl;

    /// store single track files
	std::string out_track_files = io::create_directory(OUTPUTDIR, "track_files/");

	for(auto & tracks_class_singlecell : tracks_class_singlecells){
		tracks_class_singlecell.printToFile(out_track_files);
	}

	for(auto & i : tracks_class_cluster){
        i.printToFile(out_track_files);
	}

    std::cout << " [8.3] stored all single tracks" << std::endl;

    /// store entire track files
    std::string file_tracks = io::create_directory(OUTPUTDIR, "track_files_all_in_one" );

    outputs::printDataSet(file_tracks+"/tracks.txt", tracks_class_singlecells);
    std::cout << " [8.4] stored printDataSet" << std::endl;

    /// store entire region files for each track
    std::string file_regions = io::create_directory(OUTPUTDIR, "region_files");

    outputs::printDataSetRegions(file_regions, tracks_class_singlecells);
    std::cout << " [8.5] stored printDataSetRegions" << std::endl;

    /// print total calculation time
    auto end = system_clock::now();
    const double elapsed_seconds = duration<double>(end-start).count();
    std::cout << "\n [10] End Tracking, took\t" << elapsed_seconds << " seconds" << std::endl;

    return 0;

}
