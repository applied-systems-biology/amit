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
#include <algorithm>
#include <unistd.h>
#include "IOputs.h"
#include "state_classification.h"

using std::chrono::system_clock;
using std::chrono::duration;
using json = nlohmann::json;


auto start = system_clock::now();


/************************************************************/
///                                                          /
///                           main                           /
///                                                          /
/************************************************************/
int main(int argc, char * argv[]) {

    /******************************************************************************/
    //            parse input arguments of command line and get JSON               /
    /******************************************************************************/
    json J = io::get_JSON(argc, argv);

    /// print information, like date, working directory, start time measurement
	time_t rawtime;
	struct tm *timeinfo;
	time (&rawtime);
	timeinfo = localtime(&rawtime);

    std::cout << " [0.0] Start State-Classification at Date:\t\t" << asctime(timeinfo) << std::endl;
    std::cout << " [0.1] Working directory:\t\t" << get_current_dir_name() << std::endl;

    /******************************************************************************/
    ///          assign parameter for input, output etc.                           /
    /******************************************************************************/

    std::string INPUTDIR_GRAY = J["input_gray"].get<std::string>();
    std::string INPUTDIR_BINARY = J["input_binary"].get<std::string>();
    std::string INPUTDIR_FUNGI = J["input_fungal"].get<std::string>();
    std::string INPUTDIR_DEAD = J["input_dead"].get<std::string>();
    std::string OUTPUTDIR = J["output"].get<std::string>();

    const bool FLAG_DEBUG = J["debug"].get<bool>();
    size_t N_THREADS = J["n_threads"].get<int>();

    /// print parameters to screen
    std::cout << "\t\tinput-binary directory:\t\t" << INPUTDIR_BINARY << std::endl;
    std::cout << "\t\tinput-fungiCells directory:\t" << INPUTDIR_FUNGI << std::endl;
    std::cout << "\t\tinput-deadCells directory:\t" << INPUTDIR_DEAD << std::endl;
    std::cout << "\t\toutput directory:\t\t" << OUTPUTDIR << std::endl;
    std::cout << "\t\tnumber of threads:\t\t" << N_THREADS << std::endl;

    (void) io::create_directory(OUTPUTDIR);

    /******************************************************************************/
    ///                read input images and create output directories             /
    /******************************************************************************/

    std::vector<cv::Mat> images_gray, images_binary, images_fungi, images_dead;

    io::read_images(INPUTDIR_BINARY, images_gray, cv::IMREAD_GRAYSCALE);
    io::read_images(INPUTDIR_BINARY, images_binary, cv::IMREAD_GRAYSCALE);
    io::read_images(INPUTDIR_FUNGI, images_fungi, cv::IMREAD_GRAYSCALE);
    io::read_images(INPUTDIR_DEAD, images_dead, cv::IMREAD_GRAYSCALE);

    std::cout << " [0.2] Number of loaded -brightfield- images:\t\t" << images_gray.size() << std::endl;
    std::cout << " [0.3] Number of loaded -brightfield_binary- images:\t" << images_binary.size() << std::endl;
    std::cout << " [0.4] Number of loaded -fungi_binary- images:\t\t" << images_fungi.size() << std::endl;
    std::cout << " [0.5] Number of loaded -dead_binary- images:\t\t" << images_dead.size() << std::endl;

    if (images_gray.size() == 0 || images_binary.size() == 0 || images_fungi.size() == 0 || images_dead.size() == 0){
        std::cout << "\n [INFO] One of the required images are not available. Check your directories path in the config file!" << std::endl;
        return -1;
    }

    /******************************************************************************/
    ///                     Perform State-Classification                           /
    /******************************************************************************/

    std::cout << " [2.0] Perform State-Classification" << std::endl;

    state_classification::stateClassification(images_binary, images_fungi, images_dead, OUTPUTDIR, N_THREADS, FLAG_DEBUG);
    
    /// print total calculation time
    auto end = system_clock::now();
    const double elapsed_seconds = duration<double>(end-start).count();
    std::cout << "\n [10] End StateClassification, took\t" << elapsed_seconds << " seconds" << std::endl;

    return 0;

}
