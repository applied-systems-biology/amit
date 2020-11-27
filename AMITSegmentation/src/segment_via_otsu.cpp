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

#include "segment_via_otsu.h"
#include <omp.h>
#include "IOputs.h"
#include "segment_via_ibp.h"
#include "ImageProcessingToolbox.h"

namespace IPT = ImageProcessingToolbox;


namespace otsu
{

    /**
     * Performs Otsu's Thresholding on a fluorescence image (containing fungi or dead cells).
     * Some preprocessing is done beforehand.
     *
     * Steps:
     *  - median blurring of the image using kernel size 5x5 to reduce noise
     *  - gaussian blur of the image using kernel size 5x5 and gamma 2 to smooth the image before thresholding
     *  - Otsu's thresholding
     *
     * @param src Input grey level fluorescence image
     * @return img_otsu = binary image containing thresholded cells
     */
    cv::Mat process_otsu(const cv::Mat &src, const int &MIN_REGION_SIZE, const bool &FLAG_DEBUG){

        cv::Mat img_median, img_gauss, img_otsu;
        
        // preform median filter
        cv::medianBlur(src, img_median, 5);

        // perform Gaussian filtering
        cv::GaussianBlur(img_median, img_gauss, cv::Size(5,5),2);
	
        // perform Otsu-thresholding
        cv::threshold(img_gauss, img_otsu, cv::THRESH_OTSU, 255, cv::THRESH_BINARY);

        // remove small objects, specified by the min number of pixels in an object
        cv::Mat img_rm_small_objects;
        IPT::bwareaopen(img_otsu, img_rm_small_objects, MIN_REGION_SIZE);
    
        return img_otsu;
    }

    /**
     * Calls the process_otsu() method in a multi-threaded manner.
     * The images are split evenly over the number of threads.
     * The images are chosen according to the thread's id tid.
     * The result images are saved in the output directory.
     * The percental progress is printed onto the screen.
     *
     * @param tid thread id
     */
    void otsu_segmentation(const std::vector<cv::Mat> &images, const int &N_THREADS, const int &MIN_REGION_SIZE, const std::string &OUT_BINARY, const bool &FLAG_DEBUG){

        // create directory to save control images for FLAG_DEBUG mode
        std::string OUTPUTDIR_CNT = OUT_BINARY + "/../red_binary_control/";
        if(FLAG_DEBUG)
            (void) io::create_directory(OUTPUTDIR_CNT);

        // multi-threaded processing of brightfield images
        #pragma omp parallel for firstprivate(images) num_threads(N_THREADS)
        for (size_t i = 0; i < images.size(); i++)
        {
            cv::Mat binary = otsu::process_otsu(images[i], FLAG_DEBUG, MIN_REGION_SIZE);

            if (FLAG_DEBUG) {
                cv::Mat img_cnt;
                // for debugging purpose: draw contours in source-image 
                ibp::myfindContours(binary, images[i], img_cnt);        
                // save control image
                io::write_image(OUTPUTDIR_CNT, img_cnt, i, true, std::nullopt, "jpg");
            } 

            // save binary image including PMNs
            io::write_image(OUT_BINARY, binary, i, true);
                        
        } 

    }
    
} // otsu