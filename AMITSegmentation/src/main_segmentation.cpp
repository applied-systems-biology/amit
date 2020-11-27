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
#include <unistd.h>
#include <vector>
#include <chrono>
#include <algorithm>
#include "IOputs.h"
#include "segment_via_GMM.h"
#include "segment_via_canny.h"
#include "segment_via_ibp.h"
#include "segment_via_otsu.h"
#include "segment_via_gmmFCs.h"
#include "segment_via_singleCells.h"

using std::chrono::system_clock;
using std::chrono::duration;
using json = nlohmann::json;


auto start = system_clock::now();


/************************************************************/
//                                                          //
//                           main                           //
//                                                          //
/************************************************************/
int main(int argc, char* argv[]) {
    
    /******************************************************************************/
    //            parse input arguments of command line and get JSON               /
    /******************************************************************************/
    json J = io::get_JSON(argc, argv);

    // print information, like date, working directory, start time measurement
    time_t rawtime;
  	struct tm *timeinfo;
  	time (&rawtime);
  	timeinfo = localtime(&rawtime);    

    std::cout << " [0.0] Start AMITSegmentation at Date: " << asctime(timeinfo) << std::endl;
    std::cout << " [0.1] Working directory:\t" << get_current_dir_name() << std::endl;
    
    /******************************************************************************/
    ///          assign parameter, distinguish input directory depends             /
    ///                  method and its corresponding cell-type                    /
    /******************************************************************************/
    std::string selected_method;
    
    const char *minit[] = { "gmm", "faf", "gmmFCs", "canny", "singlecells", "otsu" };
    std::vector<std::string> methods(minit, minit+6);
    
    // check for a valid method 
    if(std::find(methods.begin(), methods.end(), J["method"].get<std::string>()) != methods.end()){
        selected_method = J["method"].get<std::string>();
    } else {
        std::cout << " [ERROR] Please specify valid method: {gmm, faf, gmmFCs, canny, singlecells, otsu}" << std::endl;
        return -1;
    }

    size_t N_THREADS = J["n_threads"].get<int>();
    bool FLAG_DEBUG = J["debug"].get<bool>();

    // get the input directory
    std::string INPUTDIR = J["input"].get<std::string>();
    std::string INPUTDIR_FUNGI = J["input_gmmFCs"].get<std::string>();
        
    // create output directory
    std::string OUTPUTDIR = J["output"].get<std::string>();
    (void) io::create_directory(OUTPUTDIR);

    if(FLAG_DEBUG) 
        std::cout << " [0.3.0] Selected parameter in JSON" << J << std::endl;

    // print parameters to screen
    std::cout << "\t\tinput directory: " << INPUTDIR << std::endl;
	std::cout << "\t\toutput directory: " << OUTPUTDIR << std::endl;
    std::cout << "\t\tsegmentation-method: " << selected_method << std::endl;
	std::cout << "\t\tnumber of threads: " << N_THREADS << std::endl;
    
    /******************************************************************************/
    //read input images and create output direcorys dependently on selected method /
    /******************************************************************************/

    std::vector<cv::Mat> images, images_fungi;

    if( selected_method == "gmm" || selected_method == "faf" || selected_method == "otsu" ){

		io::read_images(INPUTDIR, images, cv::IMREAD_GRAYSCALE);
        
    } 
    else if(selected_method == "gmmFCs") {   

	    io::read_images(INPUTDIR, images, cv::IMREAD_GRAYSCALE);
        io::read_images(INPUTDIR_FUNGI, images_fungi, cv::IMREAD_COLOR);
        
        std::cout << " [0.2.1] Number of loaded brightfield-images: " << images.size() << std::endl;
        std::cout << " [0.2.2] Number of loaded green-images: " << images_fungi.size() << std::endl;

	}
    else if( selected_method == "canny" ) {

        const bool FLAG_READRGBIMAGES = false; 

		if (FLAG_READRGBIMAGES)
            io::read_images(INPUTDIR, images, cv::IMREAD_COLOR);
		else
            io::read_images(INPUTDIR, images, cv::IMREAD_GRAYSCALE);
			
    }
    else if( selected_method == "singlecells" ) {
 
        io::read_images(INPUTDIR, images, cv::IMREAD_COLOR);
			
    }
    
    std::cout << " [0.2] Number of loaded images: " << images.size() << std::endl;

    /// abort if no images loaded
    if (images.size() == 0){
        std::cout << " [WARNING] No images found: " << images.size() << std::endl;
        return -1;
    }

    /************************************************************/
    //         distinguish according to selected method          /
    /************************************************************/
    std::cout << " [0.3] Start segmentation method: " << selected_method << std::endl;
    
    if(selected_method == "gmm") { 
        
        const int C = 3;
        int N_TEMPVAR = J["number_tempvar"].get<int>();  
        int N_CLOSINGS = J["number_closing"].get<int>();
        bool CLEAR_BORDER = J["clear_border"].get<bool>();
              
        gmm::gmm_segmentation(images, N_THREADS, N_TEMPVAR, C, N_CLOSINGS, OUTPUTDIR, CLEAR_BORDER, FLAG_DEBUG);

    }

    /************************************************************************************************************/

    else if(selected_method == "faf") {
        
        std::string REMOVE_GRID = J["remove_grid"].get<std::string>();
        std::string FFT_MASK_PATH = J["fft_mask_path"].get<std::string>();
        std::cout << " [0.3.1] selected method to remove the grid: <" << REMOVE_GRID << ">" << std::endl;

        int MIN_REGION_SIZE = J["min_region_size"].get<int>();
        int SD1_KERNEL_SIZE = J["sd1_kernelSize"].get<int>();
        int SD2_KERNEL_SIZE = J["sd2_kernelSize"].get<int>();
        /// threshold should be in range of 0-255 for specified value, or <= 0 to use the mean of the corner pixel intensity
        double THRESHOLD = J["threshold_binary"].get<double>();
        std::cout << " [0.3.2] threshold for binarization:\t" << THRESHOLD << std::endl;

        int ERODE_KERNEL_SIZE = J["erode_kernelSize"].get<int>();

        /// parameter to clear all segments which are connected to the border
        bool CLEAR_BORDER = J["clear_border"].get<bool>();
        ibp::ibp_segmentation(images, N_THREADS, REMOVE_GRID, FFT_MASK_PATH, CLEAR_BORDER, OUTPUTDIR, MIN_REGION_SIZE, SD1_KERNEL_SIZE, SD2_KERNEL_SIZE,THRESHOLD, ERODE_KERNEL_SIZE, FLAG_DEBUG);
    }

    /************************************************************************************************************/

	else if(selected_method == "gmmFCs") {   

        const int C = 3;
        int N_TEMPVAR = J["number_tempvar"].get<int>();
        std::string INPUT_GMMFCs = J["input_gmmFCs"].get<std::string>();  
        bool CLEAR_BORDER = J["clear_border"].get<bool>();
        
        gmmFCs::gmmFCs_segmentation(C, images, images_fungi, N_THREADS, N_TEMPVAR, INPUT_GMMFCs, OUTPUTDIR, CLEAR_BORDER, FLAG_DEBUG);

	}

    /************************************************************************************************************/
    
    else if(selected_method == "canny") {
        
        int MORPH_OPEN_CANNY = J["morph_open"].get<int>();
        int MORPH_CLOSE_CANNY = J["morph_close"].get<int>();
        bool FLAG_MINFILTER = J["min_filter"].get<bool>();
        bool FLAG_MINFILTERPLUSMEDIAN = J["min_med_filter"].get<bool>();   
        int MIN_REGION_SIZE = J["min_region_size"].get<int>(); 
    	    
        canny::canny_segmentation(images, N_THREADS, FLAG_MINFILTER, FLAG_MINFILTERPLUSMEDIAN, MORPH_OPEN_CANNY, MORPH_CLOSE_CANNY, MIN_REGION_SIZE, OUTPUTDIR, FLAG_DEBUG);
		
	}

    /************************************************************************************************************/
    
    else if(selected_method == "otsu") {

        int MIN_REGION_SIZE = J["min_region_size"].get<int>(); 

        otsu::otsu_segmentation(images, N_THREADS, MIN_REGION_SIZE, OUTPUTDIR, FLAG_DEBUG);

    }

    /************************************************************************************************************/
    
    else if(selected_method == "singlecells"){

        singleCells::singleCells_segmentation(images, N_THREADS, OUTPUTDIR, FLAG_DEBUG);
        
    }

    // print total calculation time    
    auto end = system_clock::now();
    const double elapsed_seconds = duration<double>(end-start).count();
    std::cout << " [10] End Segment(FlatCells), took\t" << elapsed_seconds << " seconds" << std::endl;

    return 0;
}
