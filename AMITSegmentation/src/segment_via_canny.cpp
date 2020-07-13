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

#include "segment_via_canny.h"
#include <omp.h>
#include <opencv2/highgui.hpp>
#include "IOputs.h"
#include "ImageProcessingToolbox.h"
#include "segment_via_ibp.h"

namespace IPT = ImageProcessingToolbox;


namespace canny
{

    /**
     * Function corrects the contrast on an image using the passed gamma value.
     *
     * @param src image that should be contrast corrected
     * @param dst result image after gamma correction
     * @param gamma gamma value
     *
     */
    void gamma_correction(const cv::Mat &src, cv::Mat &dst, const double gamma){
        // corrects contrast on image using gamma value
        dst = cv::Mat::zeros( src.size() , src.type() );

        cv::Mat tmp, tmp2;

        src.convertTo(tmp, CV_32F, 1/255.0);
        cv::pow(tmp, gamma, tmp2);
        tmp2.convertTo(dst, CV_8U, 255);
    }
 
    /**
     * Uses Canny edge detection to find fungal cells in fluorescence images. Several steps are involved.
     *
     * Steps:
     * - Median filtering with kernel size 3x3 to reduce noise
     * - gamma correction with gamma=0.3 to increase contrast.
     * - morphological opening with kernel of size MORPH_OPEN_CANNY (default is 9) to eliminate noise particles smaller than cells. Should be adjusted to the cell size.
     * - normalization of the image to [0,255]
     * - Canny edge detection
     * - morphological closing with kernel of size MORPH_CLOSE_CANNY (default is 13) to close holes in the cell boundaries.
     * - holes in the middle of the cell clumps are filled
     *
     * @param src Input image (grey level fluorescence image)
     * @return img_canny_fill = Binary image containing the thresholded cells.
     */
    cv::Mat compute_canny_image(const cv::Mat &src, const bool &FLAG_MINFILTER, const bool &FLAG_MINFILTERPLUSMEDIAN, const int &MORPH_OPEN_CANNY, const int &MORPH_CLOSE_CANNY, const int &MIN_REGION_SIZE, const bool &FLAG_DEBUG){
        
        cv::Mat img_min, img_med, img_med_gamma, img_open;
        cv::Mat kernel;
       
        // minimum filtering is used additionally to median, if image is very noisy
        // in original version, it was used instead of median?
        if(FLAG_MINFILTER){            
            
            kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(3,3), cv::Point(1,1));
            cv::erode(src, img_min, kernel, cv::Point(-1,-1), 1);
                        
            img_min.copyTo(img_med);
        }
        else if (FLAG_MINFILTERPLUSMEDIAN){
            
            kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(3,3), cv::Point(1,1));
            cv::erode(src, img_min, kernel, cv::Point(-1,-1), 1);
            
            // perform median-filter
            cv::medianBlur(img_min, img_med, 3);
        }
        else{
            // perform only median-filter
            cv::medianBlur(src, img_med, 3);
        }

        canny::gamma_correction(img_med, img_med_gamma, 0.3);
        
        // perform morphology-opening-operation
        kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(MORPH_OPEN_CANNY, MORPH_OPEN_CANNY));
        cv::morphologyEx(img_med_gamma, img_open, cv::MORPH_OPEN, kernel);
        // std::cout << " [1.2] Used opening with ellipse of size: " << MORPH_OPEN_CANNY << std::endl;

        // normalize
        cv::normalize(img_open, img_open, 0,255, cv::NORM_MINMAX, CV_8U);
        
        cv::Mat img_canny, img_canny2, img_canny_close, img_canny_fill;

        // perform Canny - edge detector
        cv::Canny(img_open, img_canny, 0, 100, 3);
        img_canny.convertTo(img_canny2, CV_8U);

        kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(MORPH_CLOSE_CANNY, MORPH_CLOSE_CANNY));
        cv::morphologyEx(img_canny2, img_canny_close, cv::MORPH_CLOSE, kernel);
    
        // fill holes
        cv::Rect rect;
        int x = 0;
        int y = 0;

        while (rect.area() != src.rows*src.cols and (x < src.rows and y < src.cols)){
            img_canny_close.copyTo(img_canny_fill);
            cv::floodFill(img_canny_fill, cv::Point(x,y), 255, &rect, 0,0,4);
            x += 10;
        }

        cv::bitwise_not(img_canny_fill, img_canny_fill);
        cv::add(img_canny_close, img_canny_fill, img_canny_fill);

        // remove small objects, specified by the min number of pixels in an object
        cv::Mat img_rm_small_objects;
        IPT::bwareaopen(img_canny_fill, img_rm_small_objects, MIN_REGION_SIZE);

        return img_rm_small_objects; 
    }

    /**
     * Calls the canny_segmentation() method in a multi-threaded manner.
     * The images are split evenly over the number of threads.
     * The images are chosen according to the thread's id tid.
     * The result images are saved in the output directory.
     * The percental progress is printed onto the screen.
     *
     * @param tid Thread id
     * @param images read color/gray-images
     */
    void canny_segmentation(const std::vector<cv::Mat> &images, const int &N_THREADS, const bool &FLAG_MINFILTER, const bool &FLAG_MINFILTERPLUSMEDIAN, const int &MORPH_OPEN_CANNY, const int &MORPH_CLOSE_CANNY, const int &MIN_REGION_SIZE, const std::string &OUT_GREEN, const bool &FLAG_DEBUG) {
        
        // create directory to save control images for FLAG_DEBUG mode
        std::string OUTPUTDIR_CNT = OUT_GREEN + "/../green_binary_control/";
        if(FLAG_DEBUG)
            (void) io::create_directory(OUTPUTDIR_CNT);

        if(FLAG_MINFILTER){ 
            std::cout << " [1.0] Used minimum filtering" << std::endl;
        }
        else if (FLAG_MINFILTERPLUSMEDIAN) {
            std::cout << " [1.0] Used minimum and median filtering" << std::endl;
        }
        else {
            std::cout << " [1.0] Used median filtering with 3x3 window" << std::endl;
        }        
        
        std::cout << " [1.1] Minimum region size:\t" << MIN_REGION_SIZE << std::endl;

        // multi-threaded processing of fluorescence images
        #pragma omp parallel for firstprivate(images) num_threads(N_THREADS)
        for (size_t i = 0; i < images.size(); ++i)
        {
            cv::Mat binary = canny::compute_canny_image(images[i], FLAG_MINFILTER, FLAG_MINFILTERPLUSMEDIAN, MORPH_OPEN_CANNY, MORPH_CLOSE_CANNY, MIN_REGION_SIZE, FLAG_DEBUG);
                
            if (FLAG_DEBUG) {
                cv::Mat img_cnt;
                // for debugging purpose: draw contours in source-image 
                ibp::myfindContours(binary, images[i], img_cnt);        
                // save control image
                io::write_image(OUTPUTDIR_CNT, img_cnt, i, true);
            }
            
            // save binary image 
            io::write_image(OUT_GREEN, binary, i, true);
        }
    }

}
