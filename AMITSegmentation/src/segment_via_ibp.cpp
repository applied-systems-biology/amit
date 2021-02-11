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


#include "segment_via_ibp.h"
#include <cstdlib>
#include <cassert>
#include <optional>
#include <algorithm>
#include <iomanip>
#include "misaxx_convolve.h"
#include "IOputs.h"
#include "ImageProcessingToolbox.h"
#include "visualize.h"

namespace IPT = ImageProcessingToolbox;


namespace ibp
{

    /**
     * Calculate the slope for two given points with their corresponding coordinates.
     *
     * @param x0 x-coordinate of first point.
     * @param y0 y-coordinate of first point.
     * @param x1 x-coordinate of second point.
     * @param y1 y-coordinate of second point.
     * @return slope (0 means no change).
     */
    double Slope(int x0, int y0, int x1, int y1){
        return (double)(y1-y0) / (x1-x0);
    }

    /**
     * Draw a line through 2 points with the respective slope until the line connects with the image border.
     *
     * @param img image where the line will be drawn.
     * @param a first point where the line goes through.
     * @param b second point where the line goes through.
     * @param color color of the drawn line.
     */
    void fullLine(cv::Mat *img, const cv::Point& a, const cv::Point& b, const cv::Scalar& color){
        
        double slope = Slope(a.x, a.y, b.x, b.y);

        cv::Point p(0,0), q(img->cols,img->rows);

        p.y = -(a.x - p.x) * slope + a.y;
        q.y = -(b.x - q.x) * slope + b.y;

        line(*img, p, q, color, 5, 8, 0);
    }

    /**
     * Detect possible grid lines within an image with a hough transformation.
     *
     * @param src image where the grid lines supposed to be in.
     * @param dst resulting image with a binary image, where the foreground is defined by the detected lines.
     * @param minIntersection the minimum number of points that can form a line. Lines with less than this number of points are disregarded.
     * @param minLineLength minimum line length. Line segments shorter than that are rejected.
     * @param maxLineGap maximum allowed gap between points on the same line to link them.
     * @param use_canny perform canny on gray-scale image before hough transformation or not.
     * @param draw_full_line draw only the contour of the lines (draw_full line=false) or draw until the line connects with the image boundary.
     */
    void houghTransformation(const cv::Mat &src, cv::Mat &dst, const int &minIntersection, const int &minLineLength, const int &maxLineGap, const bool &use_canny,  const bool &draw_full_line) {

        cv::Mat temp;

        // check for empty dst-matrix
        if(dst.empty())
            dst = cv::Mat::zeros( src.size(), CV_32FC1 );
        
        // use Canny edge detection if choosen
        if(use_canny){
            const int kSize = 3;
            const int lowTh = 12;
            const int highTh = 12*kSize; 
            cv::Canny(src, temp, lowTh, highTh, kSize);
        }
        else {
            src.copyTo(temp);
        }    

        // Probabilistic Line Transform            
        std::vector<cv::Vec4i> lines;
        cv::HoughLinesP(temp, lines, 1, CV_PI/180, minIntersection, minLineLength, maxLineGap ); 
        // Standard Hough Transformation
        // cv::HoughLines(temp, lines, 1, CV_PI/180, minIntersection );

        for(auto l : lines) {
            
            // draw the line until it connects with the image boundary
            if(draw_full_line){
                fullLine( &dst, cv::Point(l[0], l[1]), cv::Point(l[2], l[3]), cv::Scalar(255) ); 
            }
            // draw the exact line contour
            else {
                cv::line( dst, cv::Point(l[0], l[1]), cv::Point(l[2], l[3]), cv::Scalar(255), 7, cv::LINE_AA);
            }
        
        }            
        
        dst.convertTo( dst, CV_8UC1 );
    }

    /**
     * Find contours and draw these in the input image.
     *
     * @param src_binary segmented and binary input image.
     * @param src_original original color/gray-scaled image on which the contours are drawn.
     * @param dst resulting image with a color/gray-scaled image, where contours are shown.
     * @param color color of the contours in the resulting image.
     */
    void myfindContours(const cv::Mat &src_binary, const cv::Mat &src_original, cv::Mat &dst, const cv::Scalar color) {
        
        std::vector<std::vector<cv::Point> > contours;
        std::vector<cv::Vec4i> hierarchy;

        cv::findContours( src_binary, contours, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );

        src_original.copyTo(dst);
        cv::cvtColor(dst, dst, cv::COLOR_GRAY2RGB);

        for( size_t i = 0; i< contours.size(); ++i ) 
            cv::drawContours( dst, contours, i, color, 2, 8, hierarchy, 0, cv::Point() );
        
    }

    /**
     * Perform a second hough transformation to ensure just keep the lines
     * 
     * (1) pre-processing by removoval of artifacts and objects smaller than a specified area size 
     * (2) close the remaining non-zero pixels
     * (3) perform the second hough and draw all detected lines until they connect with the image boundary
     * (4) closing operation to combine the detected lines and to prevent holes
     *
     * @param src image where the grid lines should be enhanced.
     * @param dst resulting image with a binary image, where the foreground is defined by the detected lines over the whole image.
     * @param min_area_size remove binary objects smaller than this specified area size. 
     * @param minIntersection the minimum number of points that can form a line. Lines with less than this number of points are disregarded.
     * @param minLineLength minimum line length. Line segments shorter than that are rejected.
     * @param maxLineGap maximum allowed gap between points on the same line to link them.
     */
    void enhance_gridlines(const cv::Mat &src, cv::Mat &dst, const int &min_area_size, const int &minIntersection, const int &minLineLength, const int &maxLineGap){

        cv::Mat img_bwareaopen, img_closing, img_hough;

        cv::Mat kernel = cv::getStructuringElement( cv::MORPH_CROSS, cv::Size(5,5) );

        // (1) pre-processing by removal of artifacts and objects smaller than a specified area size 
        IPT::bwareaopen(src, img_bwareaopen, min_area_size);
        
        // (2) remove the remaining non-zero pixels by morphological closing
        cv::morphologyEx(img_bwareaopen, img_closing, cv::MORPH_CLOSE, kernel);

        // (3) perform the second hough and draw all detected lines until they connect with the image boundary (CAUTION: NOT fully drawn line)
        ibp::houghTransformation(img_closing, img_hough, minIntersection, minLineLength, maxLineGap, false, false);

        // (4) closing operation to combine the detected lines and to prevent holes
        cv::morphologyEx(img_hough, dst, cv::MORPH_CLOSE, kernel, cv::Point(-1,-1), 2);

    }

    /**
     * Approximate the detection of grid lines into an image by performing a hough transformation to detect the grid lines in the image.
     *      
     * (1) specify the parameter dependent on the image size
     * (2) iterate over all parameter permutations:
     *     - start with parameter for a none sensitive (rough) detection and 
     *     - change the parameter with a nested loop towards a more precise line detection
     * (3) calculate the hough transformation for the current parameter setting
     * (4) perform a second hough transformation on the already detected lines
     * (5) calculate the number of non-zero pixels as a metric for the approximation and compare it with the previous iteration
     * (6) chose the setting with the minimum change dependent on the sensitive hough detection
     * 
     * @param src image where the grid lines should be detected.
     * @param dst resulting image with a binary image, where the foreground is defined by the detected lines over the whole image.
     * @param FLAG_DEBUG verbose parameter. 
     */
    void approx_houghLines(const cv::Mat &src, cv::Mat &dst, const bool &FLAG_DEBUG) {

        // (1) specify the parameter dependent on the image size
        std::vector<int> minIntersections( 1 );
        std::vector<int> maxLineGap( 4 );

        int minLineLength;
        int min_area_size = 1;

        // minLineLength is fixed for the size
        if (src.cols < 1500){
            // image size ~ 1024x1024
            minLineLength = 200;

            // step-parameter
            minIntersections = { 60 };
            maxLineGap = { 5, 10, 15, 20 };

            min_area_size = 1250;
        } else {
            // image size ~ 2048x2048
            minLineLength = 300;
            
            // step-parameter
            minIntersections = { 120 };
            maxLineGap = { 5, 10, 15, 20 };

            min_area_size = 2500;
        }

        // save mask with minimum change of non-zero pixels
        int predecessor_cnt = 0;
        int min_change = std::numeric_limits<int>::max();

        // (2) iterate over all parameter permutations
        for (int minInt : minIntersections) {

            for (int minLGap : maxLineGap) {

                // (3) calculate the hough transformation for the current parameter setting (no canny because input is already a binary image)
                cv::Mat img_hough;
                ibp::houghTransformation(src, img_hough, minInt, minLineLength, minLGap, false, false);

                // (4) perform a second hough transformation on the already detected lines
                cv::Mat img_hough_2nd;
                ibp::enhance_gridlines(img_hough, img_hough_2nd, min_area_size, minInt*3, minLineLength, minLGap*3);

                // (5)  calculate the number of non-zero pixels as a metric for the approximation and compare it with the previous iteration
                int current_cnt = cv::countNonZero(img_hough_2nd);

                // calculate the difference of non-zero pixels of the two images
                int absolute_diff = abs( predecessor_cnt - current_cnt );
                double relat_diff = abs( 1. - (double) current_cnt / predecessor_cnt ) * 100;

                // (6) chose the setting with the minimum change dependent on the sensitive hough detection
                if ( absolute_diff < min_change ) {
                    min_change = absolute_diff;
                    img_hough_2nd.copyTo( dst );
                }

                if (FLAG_DEBUG) {
                    std::cout << "minIntersections: " << minInt << " , maxLineGap_steps: " << minLGap  << "\tabsolute_diff:\t" << absolute_diff << "\trelative_diff:\t" << relat_diff << "\tmin_change:\t" << min_change << std::endl;
                }

                // set the current count of non-zero pixels to the predecessor for the next iteration
                predecessor_cnt = current_cnt;

            }
        }

        // post-processing by remove of artifacts and objects smaller than a specified area size
        cv::Mat img_bwareaopen;
        IPT::bwareaopen(dst, img_bwareaopen, min_area_size);

        // post-processing by remove objects which are not connected to the image border
        cv::Mat img_clearBorder;
        IPT::imclearborder(img_bwareaopen, img_clearBorder);
        cv::subtract(img_bwareaopen, img_clearBorder, dst);

    }

    /**
     * Segment a given gray-scaled based on a contrast enhancing with several additional features.
     * 
     * @param src_tmp image which should be segmented.
     * @param sd1 first special created filter to enhance the signal of cell based on their morphology.
     * @param sd2 second special created filter to enhance the signal of cell based on their morphology.
     * @param rmg method to detect possible static grid lines within the given image: {hough, fft_mask}. Every other entry means to special line detection will be used.
     * @param FFT_MASK_PATH if <rmg> = "fft_mask", this path provides the path to the directory where the corresponding mask is located.
     * @param clear_border remove objects which are connected with the image border or not.
     * @param min_region_size remove binary objects smaller than this specified area size.
     * @param OUTPUTDIR_CNT path for control images if <FLAG_DEBUG> = true.
     * @param i necessary index parameter for <OUTPUTDIR_CNT>.
     * @param FLAG_DEBUG verbose parameter. 
     * @return resulting binary image.
     */
    cv::Mat segSplitFull(const cv::Mat &src_tmp, const cv::Mat &sd1, const cv::Mat &sd2, const std::string &rmg, const std::string &FFT_MASK_PATH, const bool &clear_border, const int &min_region_size, double &mythresh, const int&erode_kernelSize, const std::string &OUTPUTDIR_CNT, const int &i, const bool &FLAG_DEBUG){

        cv::Mat src;

        ///// remove artifacts by FFT-mask
        if (rmg == "fft_mask"){

            /// read fft-mask
            cv::Mat mask = cv::imread(FFT_MASK_PATH, cv::ImreadModes::IMREAD_GRAYSCALE);
            mask.convertTo(mask, CV_32F, 1.0 / 255, 0);

            /// exception handling for images where imread failed
            cv::Mat tmp_src;
        
            if(src_tmp.empty()){
                std::cout << " fill tmp_src DONE" << std::endl;
                tmp_src = cv::Mat::zeros(mask.size(), CV_8UC1);
            } else {
                src_tmp.copyTo(tmp_src);
            }
            tmp_src.convertTo(tmp_src, CV_32F, 1.0 / 255, 0);

            if (mask.empty()) {
                std::cout << " [ERROR] segSplitFull-fft_mask: mask is empty" << std::endl;
                return cv::Mat();
            }
            
            if (FLAG_DEBUG)
                std::cout << " [INFO] fft_mask.size:\t" << mask.size() << " , path: " << FFT_MASK_PATH << std::endl;

            cv::Mat img_ifft;

            /// convolution / multiplication with mask to get rid of the grid frequencies
            misaxx_convolve fft_multiplier{};
            fft_multiplier.work(tmp_src, mask, img_ifft);

            img_ifft.convertTo(img_ifft, CV_8UC1, 255, 0);

            /// cut out artifacts lines and pad it later - later solution in fft_multiplier
            int cutout = 8;
            cv::Rect roi;
            roi.x = cutout;
            roi.y = cutout;
            roi.width = img_ifft.size().width - (cutout*2);
            roi.height = img_ifft.size().height - (cutout*2);

            cv::Mat img_cutout, tmp;
            img_cutout = img_ifft(roi);
            img_cutout.copyTo(tmp);

            cv::Mat img_pad;
            cv::copyMakeBorder(tmp, img_pad, cutout, cutout, cutout, cutout, cv::BORDER_CONSTANT, cv::Scalar(0));

            src = img_pad;
        }
        else {
            src = src_tmp;
        }

        ///// (1) image preprocessing - contrast enhancement /////
        cv::Mat Jt, J1, Jb, J;

        cv::morphologyEx(src, Jt, cv::MORPH_TOPHAT, sd2);
        cv::add(src,Jt,J1);

        cv::morphologyEx(src, Jb, cv::MORPH_BLACKHAT, sd1);
        cv::subtract(J1,Jb,J); /// CAUTION: values can be < 0

        ///// (2) noise filtration by Wiener filter /////
        cv::Mat Id;

        double noiseVariance;
        IPT::wiener2(J, Id, cv::Size(5,5), noiseVariance);

        ///// (3) standard deviation map /////
        cv::Mat Ids;

        IPT::stdfilt(Id, Ids);

        ///// (4) perform raw segmentation with optional hough transformation /////
        cv::Mat1b img_segmented;

        IPT::imbinarize(Ids, img_segmented, mythresh);

        if (FLAG_DEBUG) visual::createGUI(Ids, visual::globThreshold);

        ///// (4.1) run hough-transformation on resulting standard-deviation map /////
        cv::Mat img_raw_segmentation;

        if (rmg == "hough"){

            cv::Mat gridMask1, gridMask2, img_stdfiltering, gridArea, gridOverlay, gridMax, gridClosing, gridMedian, gridClosing2;
            cv::Mat1b gridAreaThresh;

            try {

                /// perform an approximation to detect the grid as precise as possible
                cv::Mat img_hough_approx;
                ibp::approx_houghLines(img_segmented, gridMask1, FLAG_DEBUG);

                /// invert intensities to build mask for grid area
                cv::bitwise_not(gridMask1, gridMask2);

            }
            catch(const std::exception& e) {
                std::cerr << "segSplitFull::performSecondHoughTransformation failed\n" << std::endl;
                std::cerr << e.what() << std::endl;
            }

            /// perform a filtering with a standard-deviation kernel
            IPT::stdfilter(Id, img_stdfiltering),

            /// perform overlay, caution: only info about cells, which are located over a grid line will be preserved
            cv::bitwise_and(img_stdfiltering, gridMask1, gridArea);

            /// perform overlay, caution: all info about cell, which are located over a grid line will be removed
            cv::bitwise_and(J, gridMask2, gridOverlay);

            ///// perform several operation to keep the cell-information and get rid of the noise (rest of the line)

            /// perform morphology-closing - median filtering - thresholding - morphology-closing
            cv::Mat kernel = cv::getStructuringElement( cv::MORPH_ELLIPSE, cv::Size(3,3) );
            IPT::imclose(gridArea, gridClosing, kernel, 1);

            cv::medianBlur(gridClosing, gridMedian, 3);

            double rmG_threshold = 254.;
            IPT::imbinarize(gridMedian, gridAreaThresh, rmG_threshold);

            IPT::imclose(gridAreaThresh, gridClosing2, kernel, 2);

            /// merge images with parts with and without grid area
            cv::bitwise_and(img_segmented, gridMask2, img_raw_segmentation);
            cv::bitwise_or(img_raw_segmentation, gridAreaThresh, img_raw_segmentation);

        }
        else {
            img_segmented.copyTo(img_raw_segmentation);
        }
        
        ///// (5) image-fill, morphology erosion, remove-small-objects, clear-border, majority, erosion, majority /////
        cv::Mat img_fill, img_erode, img_bwareopen, img_clearborder, img_majority, img_erode_2, img_erode_b, dst;

        /// (5.1) fill all holes in (low-contrast) cells
        IPT::imfill(img_raw_segmentation, img_fill);

        /// (5.2) morphology erosion with small rectangular shaped kernel
        cv::Mat kernelErode = cv::getStructuringElement( cv::MORPH_CROSS, cv::Size(1,1) );
        cv::erode(img_fill, img_erode, kernelErode);

        /// (5.3) remove small objects, specified by the min. number of pixels in an object
        IPT::bwareaopen(img_erode, img_bwareopen, min_region_size);

        /// (5.4) clear the objects connected to the image border if specified
        if (clear_border) {
            IPT::imclearborder(img_bwareopen, img_clearborder);    
        }
        img_bwareopen.copyTo(img_clearborder);

        /// (5.5) perform bwmorph - majority operation to set pixel to 1 if five or more pixels in its 3-by-3 neighborhood are 1s
        IPT::bwmorph(img_clearborder, img_majority, "majority");

        /// (5.6) perform a final morphology erosion with a ellipsoid-shaped kernel with specified kernel-size and subsequently area opening
        if (erode_kernelSize > 0) {
            cv::Mat kernelErode = cv::getStructuringElement( cv::MORPH_ELLIPSE, cv::Size(erode_kernelSize,erode_kernelSize) );
            cv::erode(img_majority, img_erode_2, kernelErode);
            
            IPT::bwareaopen(img_erode_2, img_erode_b, min_region_size);
        
            /// (5.7) perform again majority operation to remove generated artifacts
            IPT::bwmorph(img_erode_b, dst, "majority");    
        }
        else {
            img_majority.copyTo(dst);
        }
        
        return dst;

    }

    /**
     * Pipeline for segmentation of all images based on a contrast enhancing with parallelization possibility.
     * 
     * @param images images which should be segmented.
     * @param N_THREADS number of threads which will be used to parallelize the processing of single images.
     * @param rmg method to detect possible static grid lines within the given image: {hough, fft_mask}. Every other entry means to special line detection will be used.
     * @param FFT_MASK_PATH if <rmg> = "fft_mask", this path provides the path to the directory where the corresponding mask is located.
     * @param clear_border remove objects which are connected with the image border or not.
     * @param OUT_BINARY path where the final segmented and binary images will be stored.
     * @param min_region_size remove binary objects smaller than this specified area size.
     * @param FLAG_DEBUG verbose parameter. 
     */
    void ibp_segmentation(const std::vector<cv::Mat> &images, const int &N_THREADS, const std::string &rmg, const std::string &FFT_MASK_PATH, const bool &clear_border, const std::string &OUT_BINARY, const int &min_region_size, const int&sd1_kernelSize, const int&sd2_kernelSize, double &mythresh, const int&erode_kernelSize, const bool &FLAG_DEBUG){

        cv::Mat sd1 = IPT::create_disk(sd1_kernelSize);
        cv::Mat sd2 = IPT::create_disk(sd2_kernelSize);

        std::cout << " [0.4.0] sd1-size: " << sd1.size() << std::endl;
        std::cout << " [0.4.1] sd2-size: " << sd2.size() << std::endl;

        // create directory to save control images for FLAG_DEBUG mode
        std::string OUTPUTDIR_CNT = OUT_BINARY + "/../binary_control/";
        if(FLAG_DEBUG)
            (void) io::create_directory(OUTPUTDIR_CNT);
        
        // multi-threaded processing of brightfield images
        #pragma omp parallel for firstprivate(images) num_threads(N_THREADS)
        for(size_t i = 0; i < images.size(); ++i) {

            cv::Mat binary = ibp::segSplitFull(images[i], sd1, sd2, rmg, FFT_MASK_PATH, clear_border, min_region_size, mythresh, erode_kernelSize, OUTPUTDIR_CNT, i, FLAG_DEBUG);

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

       
} // ibp
