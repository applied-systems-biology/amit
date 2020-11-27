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

#include "visualize.h"
#include <iostream>
#include <string>
#include <vector>
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include <stdlib.h>
#include <stdio.h>


/* Function Implementations in namespace "visual" */    
namespace visual {

    /* visualization parameter */
    cv::Mat src;
    cv::Mat dst;
    int elem = 0;
    int kernel_size = 0;
    int operator_i = 0;
    int singlePrecision = 0;
    int ratio = 3;
    cv::RNG rng(12345);
    std::string window_name = "";
   

    /*
    * @function morphology_Operations
    */
    void morphology_Operations( int, void* )
    {
        // Since MORPH_X : 0,1,2,3,4,5 and 6, singlePrecision = #iterations

        cv::Mat element = cv::getStructuringElement( elem, cv::Size( 2*kernel_size + 1, 2*kernel_size+1 ), cv::Point( kernel_size, kernel_size ) );

        try {
            cv::morphologyEx( src, dst, operator_i, element, cv::Point(-1,-1), singlePrecision);
            cv::imshow( window_name, dst );
        } catch (const std::exception& e) {
            std::cout << e.what() << std::endl;
        }

    }

    /*
    * @function global thresholding
    */
    void global_Thresholding( int, void* )
    {        
        dst = src > elem;
        cv::imshow( window_name, dst );
    }

    /*
    * @function adaptive thresholding
    */
    void adaptive_Thresholding( int, void* )
    {      
        // compute ~1/8 of image size for the neighborhood, take care for odd value
        int kSizeR = 2 * std::floor(src.rows/16) + 1;
        int kSizeC = 2 * std::floor(src.cols/16) + 1;
        int kSize = std::floor((kSizeR+kSizeC)/2);       
        if (kSize % 2 == 0)
            kSize--;
        
        int offset = elem; 
        float sensitivity = ((float)singlePrecision) / 10.; // probably choose small value (<1)

            // distinguish between forground polarity is bright (true) and dark (false)
            if (operator_i == 0)
                cv::adaptiveThreshold(src, dst, offset, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY_INV, kSize, sensitivity);
            else
                cv::adaptiveThreshold(src, dst, offset, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, kSize, sensitivity);
        
        cv::imshow( window_name, dst );
    }

    /*
    * @function canny edge detection
    */
    void cannyDetector( int, void* )
    {        
        // 3. parameter: lowThreshold
        // 4. parameter: highThreshold
        // 5. parameter: kernel_size 
        
        int kSize = 3;
        cv::Canny( src, dst, elem, elem*ratio, kSize );

        cv::imshow( window_name, dst );
    }

    /*
    * @function hough transformation for detecting straight lines
    */
    void houghTransformationDetector( int, void* )
    {        
        // Declare the output variables
        cv::Mat cdst;
        
        // Edge detection
        const int kSize = 3;
        const int lowTh = 12;
        const int highTh = 12*ratio; 
        cv::Canny(src, dst, lowTh, highTh, kSize);
        
        // Copy edges to the images that will display the results in BGR
        cv::cvtColor(dst, cdst, cv::COLOR_GRAY2BGR);
        
        // perform standard Transformation operator == 1, probabilistic Transformation if operator == 0
        if (operator_i == 1) {

            // Standard Hough Line Transform
            std::vector<cv::Vec2f> lines; // will hold the results of the detection
            // operator_i = threshold: The minimum number of intersections to "*detect*" a line
            cv::HoughLines(dst, lines, 1, CV_PI/180, elem, 0, 0 ); // runs the actual detection
            
            // Draw the lines
            for(auto & i : lines) {
                float rho = i[0], theta = i[1];
                cv::Point pt1, pt2;
                double a = cos(theta), b = sin(theta);
                double x0 = a*rho, y0 = b*rho;
                pt1.x = cvRound(x0 + 1000*(-b));
                pt1.y = cvRound(y0 + 1000*(a));
                pt2.x = cvRound(x0 - 1000*(-b));
                pt2.y = cvRound(y0 - 1000*(a));
                cv::line( cdst, pt1, pt2, cv::Scalar(0,0,255), 3, cv::LINE_AA);
            }

        }
        else if (operator_i == 0) {

            // Probabilistic Line Transform
            std::vector<cv::Vec4i> linesP; // will hold the results of the detection
            // parameter 5 = minLinLength: The minimum number of points that can form a line. Lines with less than this number of points are disregarded.
            // parameter 5 = maxLineGap: The maximum gap between two points to be considered in the same line.
            cv::HoughLinesP(dst, linesP, 1, CV_PI/180, elem, singlePrecision, kernel_size); //100, 5 ); // runs the actual detection
            
            // Draw the lines
            for(auto l : linesP) {
                cv::line( cdst, cv::Point(l[0], l[1]), cv::Point(l[2], l[3]), cv::Scalar(0,0,255), 3, cv::LINE_AA);
            }

        }
        
        cv::imshow( window_name, cdst );
    }

    /*
    * @function for find contours
    */
    void findContoursFunction( int, void* )
    {
        cv::Mat canny_output;
        std::vector<std::vector<cv::Point> > contours;
        std::vector<cv::Vec4i> hierarchy;

        // detect edges using canny
        cv:Canny( src , canny_output, elem, elem*2, 3 );
        std::cout << "Canny-low threshold: " << elem << "\tCanny-high threshold: " << elem*2 << std::endl;

        // find contours
        cv::findContours( canny_output, contours, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );

        // draw contours
        dst = cv::Mat::zeros( canny_output.size(), CV_8UC3 );
        for( int i = 0; i< contours.size(); ++i ) {
            cv::Scalar color = cv::Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
            cv::drawContours( dst, contours, i, color, 2, 8, hierarchy, 0, cv::Point() );
        }

        cv::imshow( window_name, dst );
    }

    /*
    * @perform median filter
    */
    void perform_median_filter( int, void* )
    {        
        int kSize = 2*kernel_size + 1;

        cv::medianBlur(src, dst, kSize);
        
        cv::imshow( window_name, dst );
    }

    /*
     * @perform distance transformation and normalize after it
     */
    void distance_transformation( int, void* )
    {
        int distanceType = elem + 1;
        int maskType;

        // assign mask: {3, 5, Precise=0)
        if(kernel_size == 0)
            maskType = cv::DIST_MASK_PRECISE;
        if(kernel_size == 1)
            maskType = cv::DIST_MASK_3;
        else
            maskType = cv::DIST_MASK_5;

        try {
            cv::distanceTransform(src, dst, distanceType, maskType);

            /// normalize the distance image for range = {0.0, 1.0} so we can visualize
            cv::normalize(dst, dst, 0, 1.0, cv::NORM_MINMAX);

            cv::imshow( window_name, dst );
        } catch (const std::exception& e) {
            std::cout << e.what() << std::endl;
            std::cout << "kernel-size: " << kernel_size << "\tmask-type: " << maskType << std::endl;
        }

    }

    /*
     * @find moments in image
     */
    void find_moments( int, void* )
    {

        cv::Mat canny_output;
        Canny( src, canny_output, elem, elem*2, 3 );
        std::vector<std::vector<cv::Point> > contours;
        cv::findContours( canny_output, contours, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE );
        std::vector<cv::Moments> mu(contours.size() );
        for( size_t i = 0; i < contours.size(); i++ ){
            mu[i] = cv::moments( contours[i] );
        }
        std::vector<cv::Point2f> mc( contours.size() );
        for( size_t i = 0; i < contours.size(); i++ ){
            //add 1e-5 to avoid division by zero
            mc[i] = cv::Point2f( static_cast<float>(mu[i].m10 / (mu[i].m00 + 1e-5)),
                             static_cast<float>(mu[i].m01 / (mu[i].m00 + 1e-5)) );
        }

        dst = cv::Mat::zeros( canny_output.size(), CV_8UC3 );
        for( size_t i = 0; i< contours.size(); i++ ){
            cv::Scalar color = cv::Scalar( rng.uniform(0, 256), rng.uniform(0,256), rng.uniform(0,256) );
            drawContours( dst, contours, (int)i, color, 2 );
            circle( dst, mc[i], 4, color, -1 );
        }

//        std::cout << "\t Info: Area and Contour Length \n";
//        for( size_t i = 0; i < contours.size(); i++ )
//        {
//            std::cout << " * Contour[" << i << "] - Area (M_00) = " << std::fixed << std::setprecision(2) << mu[i].m00
//                 << " - Area OpenCV: " << contourArea(contours[i]) << " - Length: " << arcLength( contours[i], true ) << std::endl;
//        }

        cv::imshow( window_name, dst );
    }



    /*
    * @function createGUI for play with some morphological parameter with an UI
    */
    void createGUI(cv::Mat &img, const imgOperation &op)
    {
        src = img.clone();
        if( !src.data )
            return; 

        // distinguish between the different operations: thresholding, morphological, etc.
        if (op == imgOperation::morphology){
            // morphological parameter
            window_name = "Morphology Transformations Demo";
            int const max_operator = 7;
            int const max_elem = 2;
            int const max_kernel_size = 21;
            int const max_iterations = 10;
                    
            cv::namedWindow( window_name, cv::WINDOW_NORMAL);

            cv::createTrackbar("Operator:\n 0:Erode - 1:Dilate - 2:Opening\n 3:Closing - 4:Gradient - 5:Tophat\n 6:Blackhat - 7:Hitmiss", window_name,
                    &operator_i, max_operator, morphology_Operations );
            cv::createTrackbar( "Element:\n 0:Rect - 1:Cross - 2:Ellipse", window_name,
                &elem, max_elem, morphology_Operations );
            cv::createTrackbar( "Kernel size:\n 2n +1", window_name,
                &kernel_size, max_kernel_size, morphology_Operations );
            cv::createTrackbar( "Number of Iterations:\n", window_name,
                    &singlePrecision, max_iterations, morphology_Operations );

            morphology_Operations( 0, 0 );                       
        } 
        else if (op == imgOperation::globThreshold) {
            // threshold parameter
            window_name = "Global Thresholding";
            int const max_threshold = 255;

            cv::namedWindow( window_name, cv::WINDOW_NORMAL );

            cv::createTrackbar( "Threshold:\n", window_name, 
                &elem, max_threshold, global_Thresholding );

            global_Thresholding( 127 , 0);
        } 
        else if (op == imgOperation::adaptiveThreshold) {
            // adaptive threshold parameter
            window_name = "Adaptive Thresholding";
            int const max_offset = 255;
            int const max_sensitivity = 300; // == 300 / 10 = 30
            int const max_foreground = 1;
                    
            cv::namedWindow( window_name, cv::WINDOW_NORMAL );
            
            cv::createTrackbar( "Offset:\n", window_name,
                &elem, max_offset, adaptive_Thresholding );
            cv::createTrackbar( "Sensitivity:\n", window_name,
                &singlePrecision, max_sensitivity, adaptive_Thresholding );
            cv::createTrackbar( "Foreground:\n 0: bright - 1: dark", window_name, 
                &operator_i, max_foreground, adaptive_Thresholding );
            
            adaptive_Thresholding( 0, 0 );                       
        } 
        else if (op == imgOperation::canny) {
            // canny edge detection
            window_name = "Canny edge detection";
            int const max_threshold_low = 127;
            
            src.convertTo( src , CV_8UC1 );
            
            cv::namedWindow( window_name, cv::WINDOW_NORMAL );

            cv::createTrackbar( "Low-Threshold:\n", window_name, 
                &elem, max_threshold_low, cannyDetector );
            
            cannyDetector( 10, 0 );
        }
        else if(op == imgOperation::houghTransformation) {
            // hough transformation
            window_name = "Hough Transormation for line detection";
            int const max_operator = 1;
            int const min_num_intersections_to_detect_line = 500; 
            int const min_line_length = 1000; 
            int const max_line_gap = 200; 
            
            src.convertTo( src , CV_8UC1 );

            cv::namedWindow( window_name, cv::WINDOW_NORMAL );

            cv::createTrackbar( "Hough-Transformation-Method:\n 0: Probabilistic Line Transform \n 1: Standard Line Transform", window_name, 
                &operator_i, max_operator, houghTransformationDetector );   
            cv::createTrackbar( "Min-number of intersections-Threshold:\n", window_name, 
                &elem, min_num_intersections_to_detect_line, houghTransformationDetector ); 
            cv::createTrackbar( "Min-LineLength: The minimum number of points that can form a line: (only Probabilistic-Method)\n", window_name, 
                &singlePrecision, min_line_length, houghTransformationDetector ); 
            cv::createTrackbar( "Max-LineGap: The maximum gap between two points to be considered in the same line: (only Probabilistic-Method)\n", window_name, 
                &kernel_size, max_line_gap, houghTransformationDetector ); 

            houghTransformationDetector( 0, 0 );              
        }
        else if (op == imgOperation::findContours) {
            // find contours
            window_name = "Find Contours";
            int const max_thresh = 255;
            
            cv::namedWindow( window_name, cv::WINDOW_NORMAL );

            cv::createTrackbar( "Canny-Threshold:\n", window_name, 
                &elem, max_thresh, findContoursFunction );

            findContoursFunction( 0, 0 );
        }   
        else if (op == imgOperation::medianFiler)
        {
            // median filter
            window_name = "Median filter";
            int const max_kernel_size = 21;

            cv::namedWindow( window_name, cv::WINDOW_NORMAL );

            cv::createTrackbar( "Median-filter-kernel-size:\n 2n +1", window_name, 
                &kernel_size, max_kernel_size, perform_median_filter );

            perform_median_filter( 0, 0 );
        }
        else if (op == imgOperation::distanceTransform)
        {
            // distance transformation
            window_name = "Distance Transformation";

            int const max_elem = 7;
            int const max_kernel_size = 2;

            cv::namedWindow( window_name, cv::WINDOW_NORMAL );

            cv::createTrackbar( "Distance Metric:\n [L1, L2, C, L12, FAIR, WELSCH, HUBER]", window_name,
                                &elem, max_elem, distance_transformation );
            cv::createTrackbar( "Distance Mask:\n [ PRECISE, 3x3 ,5x5 ]", window_name,
                                &kernel_size, max_kernel_size, distance_transformation );

            distance_transformation( 0, 0 );
        }
        else if (op == imgOperation::moments)
        {
            // find moments
            window_name = "Moments";

            int const max_threshold = 255;

            cv::namedWindow( window_name, cv::WINDOW_NORMAL );

            cv::createTrackbar( "Canny-Threshold:\n", window_name,
                                &elem, max_threshold, find_moments );

            find_moments( 10, 0 );
        }
             
            
        cv::waitKey(0);
   }

    /**
     * @function show_histogram for grayscaled intensities of an image
     **/
    void show_histogram(std::string const& name, cv::Mat1b const& image) {

        // Set histogram bins count
        int bins = 256;
        int histSize[] = {bins};
        // Set ranges for histogram bins
        float lranges[] = {0, 256};
        const float* ranges[] = {lranges};
        // create matrix for histogram
        cv::Mat hist;
        int channels[] = {0};

        // create matrix for histogram visualization
        int const hist_height = 256;
        cv::Mat3b hist_image = cv::Mat3b::zeros(hist_height, bins);

        cv::calcHist(&image, 1, channels, cv::Mat(), hist, 1, histSize, ranges, true, false);

        double max_val=0;
        minMaxLoc(hist, 0, &max_val);

        // visualize each bin
        for(int b = 0; b < bins; b++) {
            float const binVal = hist.at<float>(b);
            int   const height = cvRound(binVal*hist_height/max_val);
            cv::line
                    ( hist_image
                            , cv::Point(b, hist_height-height), cv::Point(b, hist_height)
                            , cv::Scalar::all(255)
                    );
        }
        cv::imshow(name, hist_image);
    }

}




