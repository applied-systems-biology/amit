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

/*
This namespace contains several functions of a script language image processing toolbox that were transferred to C++ using the OpenCV library

The most parameter are first the source (src) image and secondly the destination (dst) image
*/


#include "ImageProcessingToolbox.h"
#include <stdlib.h>
#include <math.h> 
#include <assert.h>
#include <algorithm>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>


namespace ImageProcessingToolbox
{

    std::string type2str(int type) {
        std::string r;

        uchar depth = type & CV_MAT_DEPTH_MASK;
        uchar chans = 1 + (type >> CV_CN_SHIFT);

        switch ( depth ) {
            case CV_8U:  r = "8U"; break;
            case CV_8S:  r = "8S"; break;
            case CV_16U: r = "16U"; break;
            case CV_16S: r = "16S"; break;
            case CV_32S: r = "32S"; break;
            case CV_32F: r = "32F"; break;
            case CV_64F: r = "64F"; break;
            default:     r = "User"; break;
        }

        r += "C";
        r += (chans+'0');

        return r;
    }

    /**
    * structure elements for morphological operations
    *
    * @return kernel in format: CV_8U
    */
    cv::Mat createKernel(const std::string kernelType){
        cv::Mat kernel;
        if (kernelType == "disk3"){
            kernel = cv::Mat::ones( cv::Size(5,5), CV_32F  ); // =strel('disk', 3);
        }
        else if(kernelType == "disk5"){
            kernel = cv::Mat::ones( cv::Size(9,9), CV_32F ); // =strel('disk', 5);
            int N = kernel.rows-1;
            // top left
            kernel.at<float>(cv::Point(0,0)) = 0;
            kernel.at<float>(cv::Point(0,1)) = 0;
            kernel.at<float>(cv::Point(1,0)) = 0;
            // top right
            kernel.at<float>(cv::Point(0,N-1)) = 0;
            kernel.at<float>(cv::Point(0,N)) = 0;
            kernel.at<float>(cv::Point(1,N)) = 0;
            // bottom left
            kernel.at<float>(cv::Point(N-1,0)) = 0;
            kernel.at<float>(cv::Point(N,0)) = 0;
            kernel.at<float>(cv::Point(N,1)) = 0;
            // bottom right
            kernel.at<float>(cv::Point(N-1,N)) = 0;
            kernel.at<float>(cv::Point(N,N-1)) = 0;
            kernel.at<float>(cv::Point(N,N)) = 0;
        }
        else if(kernelType == "disk10"){
            kernel = cv::Mat::ones( cv::Size(20,20), CV_32F ); // =strel('disk', 10);
            int N = kernel.rows-1;
            // top left
            kernel.at<float>(cv::Point(0,0)) = 0; kernel.at<float>(cv::Point(0,1)) = 0; kernel.at<float>(cv::Point(0,2)) = 0; kernel.at<float>(cv::Point(0,3)) = 0;
            kernel.at<float>(cv::Point(1,0)) = 0; kernel.at<float>(cv::Point(1,1)) = 0; kernel.at<float>(cv::Point(1,2)) = 0;
            kernel.at<float>(cv::Point(2,0)) = 0; kernel.at<float>(cv::Point(2,1)) = 0;
            kernel.at<float>(cv::Point(3,0)) = 0;
            // top right
            kernel.at<float>(cv::Point(0,N-3)) = 0; kernel.at<float>(cv::Point(0,N-2)) = 0; kernel.at<float>(cv::Point(0,N-1)) = 0; kernel.at<float>(cv::Point(0,N)) = 0;
            kernel.at<float>(cv::Point(1,N-2)) = 0; kernel.at<float>(cv::Point(1,N-1)) = 0; kernel.at<float>(cv::Point(1,N)) = 0;
            kernel.at<float>(cv::Point(2,N-1)) = 0; kernel.at<float>(cv::Point(2,N)) = 0;
            kernel.at<float>(cv::Point(3,N)) = 0;
            // bottom left
            kernel.at<float>(cv::Point(N-3,0)) = 0;
            kernel.at<float>(cv::Point(N-2,0)) = 0; kernel.at<float>(cv::Point(N-2,1)) = 0;
            kernel.at<float>(cv::Point(N-1,0)) = 0; kernel.at<float>(cv::Point(N-1,1)) = 0; kernel.at<float>(cv::Point(N-1,2)) = 0;
            kernel.at<float>(cv::Point(N,0)) = 0; kernel.at<float>(cv::Point(N,1)) = 0; kernel.at<float>(cv::Point(N,2)) = 0; kernel.at<float>(cv::Point(N,3)) = 0;
            // bottom right
            kernel.at<float>(cv::Point(N-3,N)) = 0;
            kernel.at<float>(cv::Point(N-2,N)) = 0; kernel.at<float>(cv::Point(N-2,N-1)) = 0;
            kernel.at<float>(cv::Point(N-1,N)) = 0; kernel.at<float>(cv::Point(N-1,N-1)) = 0; kernel.at<float>(cv::Point(N-1,N-2)) = 0;
            kernel.at<float>(cv::Point(N,N)) = 0; kernel.at<float>(cv::Point(N,N-1)) = 0; kernel.at<float>(cv::Point(N,N-2)) = 0; kernel.at<float>(cv::Point(N,N-3)) = 0;
        }
        else if(kernelType == "line3_135"){
            kernel = cv::Mat::zeros( cv::Size(3,3) , CV_32F ); //=strel('line',3,135);
            kernel.at<float>(cv::Point(0,0)) = 1;
            kernel.at<float>(cv::Point(1,1)) = 1;
            kernel.at<float>(cv::Point(2,2)) = 1;
        }
        else if(kernelType == "line3_45") {
            kernel = cv::Mat::zeros( cv::Size(3,3) , CV_32F ); //=strel('line',3,45);
            kernel.at<float>(cv::Point(0,2)) = 1;
            kernel.at<float>(cv::Point(1,1)) = 1;
            kernel.at<float>(cv::Point(2,0)) = 1;
        }
        else if (kernelType == "diamond1") {
            kernel = cv::getStructuringElement( cv::MORPH_CROSS, cv::Size(3,3) );
        }
        else {
            std::cout << "please specify valid kernel type: [disk3, disk5, line3_135, line3_45, diamond1]" << std::endl;
            kernel = cv::Mat::ones( cv::Size(9,9), CV_8U ); // default
        }
        cv::Mat output;
        kernel.convertTo( output, CV_8U );
        return output;
    }

    bool contourTouchesImageBorder(std::vector<cv::Point> &contour, cv::Size &imageSize) {
        cv::Rect bb = cv::boundingRect(contour);

        bool retval = false;

        int xMin, xMax, yMin, yMax;

        xMin = 0;
        yMin = 0;
        xMax = imageSize.width - 1;
        yMax = imageSize.height - 1;

        // Use less/greater comparisons to potentially support contours outside of 
        // image coordinates, possible future workarounds with cv::copyMakeBorder where
        // contour coordinates may be shifted and just to be safe.
        if( bb.x <= xMin || 
            bb.y <= yMin ||
            bb.width >= xMax ||
            bb.height >= yMax)
        {
            retval = true;
        }

        return retval;
    }


    /**
     *  extract unique values in cv::Mat 
     */
    std::vector<float> unique(const cv::Mat& src, bool sort) {
        cv::Mat input;
        src.convertTo( input, CV_32F );

        if (input.channels() > 1 || input.type() != CV_32F) 
        {
            std::cerr << "unique !!! Only works with CV_32F 1-channel Mat" << std::endl;
            return std::vector<float>();
        }

        std::vector<float> out;
        for (int y = 0; y < input.rows; ++y)
        {
            const float* row_ptr = input.ptr<float>(y);
            for (int x = 0; x < input.cols; ++x)
            {
                float value = row_ptr[x];

                if ( std::find(out.begin(), out.end(), value) == out.end() )
                    out.push_back(value);
            }
        }

        if (sort)
            std::sort(out.begin(), out.end());

        return out;
    }

    /**
     *  compute gradient magnitude of an image
     */
    void imgradient(const cv::Mat &src, cv::Mat &magnitude) {
        // find x and y gradients
        cv::Mat sobelx, sobely;
        cv::Sobel(src, sobelx, CV_64F, 1, 0);
        cv::Sobel(src, sobely, CV_64F, 0, 1);

        // compute the magnitude
        cv::magnitude(sobelx, sobely, magnitude);

        // the following (preudo)-code compute the direction (angle) of the gradients of the image 
        // angle = np.arctan2(sobely, sobelx) * (180 / np.pi) // = python
        // angle = atan2(sobely, sobelx) * (180 / PI);
    }

    /**
     *  Burn the binary image into the 2-D grayscale image, 
     *  choosing the color to be used.
     */
    void imoverlay(const cv::Mat &src, const cv::Mat &mask, cv::Mat &dst, const std::string color){
        cv::Mat1b src_temp, mask_temp;
        
        // convert to 8-bit color image (3-channel)
        src.convertTo( src_temp, CV_8UC1 );
        mask.convertTo( mask_temp, CV_8UC1 );
        cv::cvtColor(src_temp, src_temp, cv::COLOR_GRAY2RGB );
        cv::cvtColor(mask_temp, mask_temp, cv::COLOR_GRAY2RGB );

        // find all white pixels (255,255,255)-value
        cv::Mat mask_white;
        cv::inRange(mask_temp, cv::Scalar(254,254,254), cv::Scalar(255,255,255), mask_white);
        
        // set all white pixels to specified color
        if (color == "blue"){
            dst.setTo( cv::Scalar(255,0,0) , mask );
        }
        else if (color == "red") {
            dst.setTo( cv::Scalar(0,0,255) , mask );
        }        
        else if (color == "yellow") {
            dst.setTo( cv::Scalar(0,255,255), mask );
        }        
        else
        {
            std::cout << "ERROR: pleasy choose one of the following colors: { blue, red, yellow }" << std::endl;
        }
        cv::bitwise_and(mask_temp, src_temp, dst);
    }

    /**
     *  Fills holes in the input binary image. In this syntax, 
     *  a hole is a set of background pixels that cannot be reached 
     *  by filling in the background from the edge of the image.
     */
    void imfill(const cv::Mat &src, cv::Mat &dst){
        
        // Threshold. Set values equal to or above 220 to 0. Set values below 220 to 255.
        // https://www.learnopencv.com/filling-holes-in-an-image-using-opencv-python-c/
        
        cv::Mat im_th, im_floodfill_inv;
        //threshold(src, im_th, 220, 255, cv::THRESH_BINARY_INV);
        
        // floodfill from point (0, 0)
        cv::Mat im_floodfill = src.clone();

        // to avoid artefacts at the border and floodfill the whole image, set first column to 0
        for (size_t i = 0; i < im_floodfill.cols; i++)        
            im_floodfill.at<u_int8_t>(cv::Point(i,0)) = 0;
        
        floodFill(im_floodfill, cv::Point(0,0), cv::Scalar(255));
                
        // invert floodfilled image
        bitwise_not(im_floodfill, im_floodfill_inv);
               
        // combine the two images to get the foreground.
        dst = (src | im_floodfill_inv);

    }

    /**
     * Suppresses structures in image I that are lighter 
     * than their surroundings and that are connected to the image border.
     * Use this function to clear the image border.
     */   
    void imclearborder(const cv::Mat &src, cv::Mat &dst, const int &conn){
        // requires a binary src image
        dst = src > 127;
        
        cv::rectangle(dst, cv::Rect(0, 0, dst.cols, dst.rows), cv::Scalar(255));
        cv::floodFill(dst, cv::Point(0, 0), cv::Scalar(0));

        dst.convertTo( dst , CV_8UC1 );

        /* // given a black and white image, first find all of its contours; CV_8UC1 - type necessary 
        cv::Mat src_temp = src.clone();
        src_temp.convertTo( src_temp , CV_8UC1 );

        std::vector<std::vector<cv::Point> > contours;
        std::vector<cv::Vec4i> hierarchy;
        cv::findContours(src_temp, contours, hierarchy, CV_RETR_LIST, cv::CHAIN_APPROX_SIMPLE);
        
        int width = src.cols;
        int height = src.rows;

        // ID list of contours that touch the border
        std::vector<int> contourList;

        // for each counter
        for(size_t i = 0; i < contours.size(); ++i)
        {
            // get the i'th contour
            std::vector<cv::Point> cnt = contours[i];

            // look at each point in the contour
            for(size_t j = 0; j < cnt.size(); ++j)
            {
                cv::Point pt = cnt[j];
                int rowCnt = pt.x;
                int colCnt = pt.y;

                // if point is within the radius of the border, remove contour
                bool check1 = (rowCnt >= 0 && rowCnt < conn) || (rowCnt >= height-1-conn && rowCnt < height);
                bool check2 = (colCnt >= 0 && colCnt < conn) || (colCnt >= width-1-conn && colCnt < width);

                if (check1 || check2)
                    contourList.push_back(i);                
            }       
        }
        
        for(size_t c = 0; c < contourList.size(); ++c)
        {            
            try
            {
                cv::drawContours(src_temp, contours, c, cv::Scalar(255), -1);
            }
            catch(const std::exception& e)
            {
                //std::cerr << e.what() << '\n';
                //std::cout << "OpenCV Error: Assertion failed (0 <= contourIdx && contourIdx < (int)last) in drawContours: " << contourList[c] << std::endl;
                continue;
            }            
        }
        // copy and convert result image to dst
        //src_temp.convertTo( src_temp , CV_32F );
        src_temp.copyTo(dst); */
    }

    /**
     * perform imerode, erodes the grayscale or binary image with a 
     * structuring element object or array.
     */
    void imerode(const cv::Mat &src, cv::Mat &dst, const cv::Mat &kernel, const int &iterations) {
        cv::morphologyEx(src, dst, cv::MORPH_ERODE, kernel, cv::Point(-1,-1), iterations, cv::BORDER_DEFAULT); 
        dst.convertTo( dst , CV_8UC1, 255, 0);
    }

    /**
     * perform imdilate, morphology operation dilate
     */
    void imdilate(const cv::Mat &src, cv::Mat &dst, const cv::Mat &kernel, const int &iterations) { 
        cv::morphologyEx(src, dst, cv::MORPH_DILATE, kernel, cv::Point(-1,-1), iterations, cv::BORDER_DEFAULT);      
        dst.convertTo( dst , CV_8UC1, 255, 0);
    } 

    /**
     * perform imclose, morphological closing on the grayscale or binary image
     * is a dilation followed by an erosion, using the same structuring element for both operations.
     */
    void imclose(const cv::Mat &src, cv::Mat &dst, const cv::Mat &kernel, const int &iterations){
        cv::morphologyEx(src, dst, cv::MORPH_CLOSE, kernel, cv::Point(-1,-1), iterations, cv::BORDER_DEFAULT);
    }

    /**
     * perform imopen, morphological opening on the grayscale or binary image
     */
    void imopen(const cv::Mat &src, cv::Mat &dst, const cv::Mat &kernel, const int &iterations){
        cv::morphologyEx(src, dst, cv::MORPH_OPEN, kernel, cv::Point(-1,-1), iterations, cv::BORDER_DEFAULT);
    }

    /**
     * Adjust image intensity values or colormap = imadjust 
     * 
     * @param src input CV_8UC1 image
     * @param dst output CV_8UC1 imge
     * @param tol tolerance, from 0 to 100.
     * @param in  src image bounds
     * @param out dst image buonds
     */
    void imadjust(const cv::Mat1b &src, cv::Mat1b &dst, int tol, cv::Vec2i in, cv::Vec2i out){
        dst = src.clone();

        tol = cv::max(0, cv::min(100, tol));

        if (tol > 0)
        {
            // Compute in and out limits
            // Histogram
            std::vector<int> hist(256, 0);
            for (size_t r = 0; r < src.rows; ++r) {
                for (size_t c = 0; c < src.cols; ++c) {
                    hist[src(r,c)]++;
                }
            }

            // Cumulative histogram
            std::vector<int> cum = hist;
            for (size_t i = 1; i < hist.size(); ++i) 
                cum[i] = cum[i - 1] + hist[i];
            
            // Compute bounds
            int total = src.rows * src.cols;
            int low_bound = total * tol / 100;
            int upp_bound = total * (100-tol) / 100;
            in[0] = distance(cum.begin(), lower_bound(cum.begin(), cum.end(), low_bound));
            in[1] = distance(cum.begin(), lower_bound(cum.begin(), cum.end(), upp_bound));

        }

        // Stretching
        float scale = float(out[1] - out[0]) / float(in[1] - in[0]);
        for (int r = 0; r < dst.rows; ++r)
        {
            for (int c = 0; c < dst.cols; ++c)
            {
                int vs = cv::max(src(r, c) - in[0], 0);
                int vd = cv::min(int(vs * scale + 0.5f) + out[0], out[1]);
                dst(r, c) = cv::saturate_cast<uchar>(vd);
            }
        }
    }

    /**
     * distinguish between global and adaptive thresholding = imbinarize 
     * 
     * @param src input CV_8UC1 image
     * @param dst output CV_8UC1 imge
     * @param globThreshold normalized threshold in range [0,1], default -1. means local thresholding will be performed
     * @param foreground indicate that the foreground is darker than the background
     */
    void imbinarize(const cv::Mat &src, cv::Mat1b &dst, const double &globThreshold, const bool &foregroundBright, const int &offset, const float &sensitivity){
        if (globThreshold > 0){
            // global thresholding
            double threshold = 255 * globThreshold;;
            dst = src > threshold;
        }
        else
        {
            // global thresholding
            // compute ~1/8 of image size for the neighborhood, take care for odd value
            int kSizeR = 2 * std::floor(src.rows/16) + 1;
            int kSizeC = 2 * std::floor(src.cols/16) + 1;
            int kSize = std::floor((kSizeR+kSizeC)/2);       
            if (kSize % 2 == 0)
                kSize--;
                    
            // distinguish between forground polarity is bright (true) and dark (false)
            if (foregroundBright)
                cv::adaptiveThreshold(src, dst, offset, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY_INV, kSize, sensitivity);
            else
                cv::adaptiveThreshold(src, dst, offset, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, kSize, sensitivity);     
        }  
    }

    /**
     * Subtract one image (subtracting) from another (minuating)
     */
    void imsubtract(const cv::Mat &minuating, const cv::Mat &subtracting, cv::Mat &dst) {
        cv::subtract(minuating, subtracting, dst);
    }

    /**
     * Add two images by their summands
     */
    void imadd(const cv::Mat &summand1, const cv::Mat &summand2, cv::Mat &dst) {
        cv::add(summand1, summand2, dst);
    }

    /**
     * Label connected components in a binary image and returns a matrix dst, 
     * of the same size as src, containing labels for the connected objects in src. 
     * Default neighborhood is 8, which specifies 8-connected objects
     */
    void bwlabel(const cv::Mat1b &src, cv::Mat &dst, int &nLabels, int connectivity){

        // get connected components     
        nLabels = cv::connectedComponents( src, dst, connectivity, CV_16U);

        /// to visualize all labels
//        cv::Mat label;
//        int n_la = connectedComponents(src, label, 8, CV_16U);
//
//        cv::Mat seeMyLabels;
//        normalize(label, seeMyLabels, 0, 255, cv::NORM_MINMAX, CV_8U);

        // std::cout << "Number of connected components/labels = " << nLabels << std::endl;
                
    }

    /**
     * Morphological operations on binary image, input AND output is CV_8UC1
     */
    void bwmorph(const cv::Mat &src, cv::Mat &dst, const std::string operation){

        cv::Mat src_temp = src.clone();
            if ( src_temp.type() != CV_8UC1 )
                src_temp.convertTo( src_temp, CV_8UC1 );

        if (operation == "hbreak") {
            // build identity kernel
            cv::Mat1b ones = cv::Mat1b::ones( cv::Size(16,1) );
            cv::Mat1b zeros = cv::Mat1b::zeros( cv::Size(16,1) );
            cv::Mat kernel_temp, kernel;
            cv::hconcat(ones, zeros, kernel_temp);
            cv::repeat(kernel_temp, 1, 16, kernel);

            // the 2 exceptions
            kernel.at<int>(cv::Point(382)) = 0; 
            kernel.at<int>(cv::Point(472)) = 0; 
            
            // apply lookup-table-kernel on image, rescale and invert intensities
            bwlookup(src_temp, dst, kernel);
            dst = dst * 255;
            dst = cv::Scalar::all(255) - dst;                
        }
        else if(operation == "spur"){
            // remove spur pixels: build identity kernel
            cv::Mat1b ones = cv::Mat1b::ones( cv::Size(16,1) );
            cv::Mat1b zeros = cv::Mat1b::zeros( cv::Size(16,1) );
            cv::Mat kernel_temp, kernel;
            cv::hconcat(ones, zeros, kernel_temp);
            cv::repeat(kernel_temp, 1, 16, kernel);

            // 4 qualifying patterns
            kernel.at<int>(cv::Point(18)) = 0; 
            kernel.at<int>(cv::Point(21)) = 0; 
            kernel.at<int>(cv::Point(81)) = 0; 
            kernel.at<int>(cv::Point(273)) = 0;

            // apply lookup-table-kernel on image
            bwlookup(src_temp, dst, kernel);
            dst = dst * 255;
            dst = cv::Scalar::all(255) - dst;            
        }
        else if (operation == "erode") {
            // 'erode' performs erosion using the structuring element ones(3) and
            cv::Mat kernel = cv::Mat::ones( cv::Size(3,3), CV_32F   );
            cv::erode(src, dst, kernel, cv::Point(-1,-1), 1, cv::BORDER_CONSTANT, cv::Scalar(0));
        }
        else if(operation == "majority") {
            // Sets a pixel to 1 if five or more pixels in its 3-by-3 neighborhood are 1s; otherwise, it sets the pixel to 0
            // first perform boxFilter, than threshold each perceptive field by x_ij >= 5 / 9
            cv::Mat dst_temp;
            cv::boxFilter(src_temp, dst_temp, -1, cv::Size(3,3) , cv::Point(-1,-1) , true, cv::BORDER_ISOLATED);
            
            dst = cv::Mat::zeros( src.size(), CV_8UC1 ); 
                        
            // global threshold at ( 5/9 )
            double threshold = 5. / 9. ;
            dst = dst_temp > threshold;

            dst.convertTo(dst, CV_8UC1 );            
        }
        else if (operation == "remove") {
            cv::Mat kernel = cv::Mat::ones( cv::Size(3,3), CV_32F );
            cv::morphologyEx(src, dst, cv::MORPH_GRADIENT, kernel, cv::Point(-1,-1), 1, cv::BORDER_DEFAULT);
        }      
        else if (operation == "dilate") {
            // 'dilate' performs dilation using the structuring element ones(3) and
            cv::Mat kernel = cv::Mat::ones( cv::Size(3,3), CV_32F );
            cv::dilate(src, dst, kernel, cv::Point(-1,-1), 1, cv::BORDER_CONSTANT, cv::Scalar(0));
        }
        else if (operation == "shrink") {
            cv::Mat src_gray = src.clone();
            // convert image to gray and blur it
            cv::cvtColor( src, src_gray, cv::COLOR_BGR2GRAY );
            cv::blur( src_gray, src_gray, cv::Size(3,3) );
        
            cv::Mat canny_output;
            std::vector<std::vector<cv::Point> > contours;
            std::vector<cv::Vec4i> hierarchy;

            // detect edges using canny
            double thresh = 100.;
            cv::Canny( src_gray, canny_output, thresh, thresh*2, 3 );
            // find contours
            cv::findContours( canny_output, contours, hierarchy, cv::RetrievalModes::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );

            // get the moments
            std::vector<cv::Moments> mu(contours.size() );
            for( size_t i = 0; i < contours.size(); ++i ) 
                mu[i] = moments( contours[i], false );             

            // get the mass centers:
            std::vector<cv::Point2f> mc( contours.size() );
            for( size_t i = 0; i < contours.size(); ++i )
                mc[i] = cv::Point2f( mu[i].m10/mu[i].m00 , mu[i].m01/mu[i].m00 ); 

            // draw contours
            cv::RNG rng(12345);
            cv::Mat drawing = cv::Mat::zeros( canny_output.size(), CV_8UC3 );
            for( int i = 0; i < contours.size(); i++ ) {
                cv::Scalar color = cv::Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
                cv::drawContours( drawing, contours, i, color, 2, 8, hierarchy, 0, cv::Point() );
                circle( drawing, mc[i], 4, color, -1, 8, 0 );
            }

            /// Show in a window
            // namedWindow( "Contours", CV_WINDOW_AUTOSIZE );
            // imshow( "Contours", drawing );

            /// Calculate the area with the moments 00 and compare with the result of the OpenCV function
            //printf("\t Info: Area and Contour Length \n");
            for( size_t i = 0; i < contours.size(); ++i ) {
                // printf(" * Contour[%d] - Area (M_00) = %.2f - Area OpenCV: %.2f - Length: %.2f \n", i, mu[i].m00, contourArea(contours[i]), arcLength( contours[i], true ) );
                cv::Scalar color = cv::Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
                cv::drawContours( drawing, contours, i, color, 2, 8, hierarchy, 0, cv::Point() );
                cv::circle( drawing, mc[i], 4, color, -1, 8, 0 );
            }        
            drawing.copyTo(dst);
        }                  
        else {
            std::cout << "ERROR: Please choose one of the following operations: { hbreak, spur, erode, majority, remove, dilate, shrink }" << std::endl;       
            dst = cv::Mat::zeros( src.size(), src.type() );
        }
    }

    /**
     * Remove small objects from binary image by removing all connected 
     * components (objects) that have fewer than p pixels.
     * This operation is known as an area opening
     * conn = 4, means: two-dimensional four-connected neighborhood
     * 
     * @param src input CV_8UC1 image
     * @param dst output CV_8UC1 imge
     * @param p minimum number of pixels
     */
    void bwareaopen(const cv::Mat &src, cv::Mat &dst,  const int &p){   

        dst = cv::Mat::zeros(src.size(), CV_8UC1);

        std::vector<std::vector<cv::Point>> contours;
        std::vector<cv::Vec4i> hierarchy;
        cv::findContours(src, contours, hierarchy, cv::RetrievalModes::RETR_LIST, cv::CHAIN_APPROX_SIMPLE);
        
        // for each contour, determine its total occupying area
        for(size_t i = 0; i < contours.size(); ++i)
        {
            double area = fabs(contourArea(cv::Mat(contours[i])));

            if( fabs( area ) >= (double) p){
                cv::drawContours(dst, contours, i, cv::Scalar(255), cv::FILLED);
            }

        }
    }

    /**
     * Measure properties of image regions. Returns measurements for the set of properties
     * specified by properties for each 8-connected component (object) in the binary image.
     */
    void regionprops(const cv::Mat &src, std::vector<std::vector<cv::Point>> &object_coordinates, const std::string &method) {
        // required input format: CV_8UC1   
        cv::Mat src_temp = src.clone();
        src_temp.convertTo( src_temp , CV_8UC1 );
                
        if (method == "PixelIdxList") {
            // variante 1: connected components (recommended)
            cv::Mat labels;
            int connectivity = 8;
            int label_count = cv::connectedComponents(src_temp, labels, connectivity, CV_32S);
            
            object_coordinates.clear();
            object_coordinates.resize(label_count-1);
            
            for (int y = 0; y < labels.rows; y++) {
                const int* row = labels.ptr<int>(y);

                for (int x = 0; x < labels.cols; x++) {
                    // skip background
                    if (row[x] > 0)
                    {
                        object_coordinates[row[x]-1].push_back( cv::Point(x,y) );
                    }                                     
                }                
            }

            // remove objects which touch the border of the image
            // for (size_t i = 0; i < object_coordinates.size(); i++)
            // {
            //     cv::Rect bounding_rect = cv::boundingRect(object_coordinates[i]);
            //     cv::Rect test_rect = bounding_rect & cv::Rect(1, 1, src_temp.cols - 2, src_temp.rows - 2);
            //     if (bounding_rect != test_rect)
            //     {
            //         cv::drawContours(src_temp, object_coordinates, (int)i, cv::Scalar(0),-1);
            //     }
            // }   

            // variante 2: find contours (not recommended) src must be binary image, just returns the contour-pixels
            // std::vector<cv::Vec4i> hierarchy;            
            // cv::findContours( src_temp, object_coordinates, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );
        }
        else
        {
            std::cout << "ERROR: Please choose one of the following methods: { PixelIdxList }" << std::endl; 
        }        
    }

    /**
     * Applies contrast-limited adaptive histogram equalization (CLAHE)
     * 
     * @param clipLimit is a contrast factor that prevents oversaturation of the image specifically in homogeneous areas
     */
    void adapthisteq(const cv::Mat &src, cv::Mat &dst, const double clipLimit, const std::string distribution){
        int gridSize = 8;

        cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE();
        clahe->setClipLimit(clipLimit);
        clahe->setTilesGridSize(cv::Size(gridSize,gridSize));
        
        clahe->apply(src, dst);
    }

    /**
     * This function performs an standard deviation filter 
     * on a specified neighborhood on normalized image [0,1]
     */
    void stdfilt(const cv::Mat &src, cv::Mat &dst){
        // define the output matrix
        dst = src.clone();
        // find out the mean image
        cv::blur(src, dst, cv::Size(3,3));
        // substract mean image
        cv::absdiff(dst, src, dst);
        // square to get std matrix
        cv::pow(dst, 2.0, dst);
    }

    void bwlookup(const cv::Mat & in, cv::Mat & out, const cv::Mat & lut, int bordertype, cv::Scalar px ) {
        if ( in.type() != CV_8UC1 )
            CV_Error( cv::Error::StsError, "er");
        if ( lut.type() != CV_8UC1 || lut.rows*lut.cols!=512 || !lut.isContinuous() )
            CV_Error( cv::Error::StsError, "lut size != 512" );
        if ( out.type() != in.type() || out.size() != in.size() ) 
            out = cv::Mat( in.size(), in.type() );        

        const unsigned char * _lut = lut.data;
        cv::Mat t;
        cv::copyMakeBorder( in,t,1,1,1,1,bordertype,px);
        const int rows=in.rows+1;
        const int cols=in.cols+1;
        for ( int y=1;y<rows;++y)
        {
            for ( int x=1;x<cols;++x)
            {
                int L = 0;
                const int jmax=y+1;
                //#if 0 // row-major order
                if(false){
                    for ( int j=y-1, k=1; j<=jmax; ++j, k<<=3 )
                    {
                        const unsigned char * p = t.ptr<unsigned char>(j) + x-1;
                        for (size_t u=0;u<3;++u )
                        {
                            if ( p[u] )
                                L += (k<<u);
                        }
                    }
                }
                //#else // column-major order
                else {
                    for ( int j=y-1, k=1; j<=jmax; ++j, k<<=1 )
                    {
                        const unsigned char * p = t.ptr<unsigned char>(j) + x-1;
                        for (size_t u=0;u<3;++u )
                        {
                            if ( p[u] )
                                L += (k<<3*u);
                        }
                //#endif
                    }
                }
                out.at<unsigned char>(y-1,x-1)=_lut[ L ];
            }
        }
    } 

    /**
     * This function applies to the src image the adaptive Wiener filter
     * and store the result in the dst image. 
     * The formula that will be apply is the following one: 
     * dst(x, y) = u + max(0, s^2 - v^2)(src(x, y) - u) / max(s^2, v^2) 
     * https://github.com/prittt/AdaptiveWienerFilter/blob/master/WienerFilter.h
     */
    void wiener2(const cv::Mat &src, cv::Mat &dst, const cv::Size &kSize, double &noiseVariance) {
        assert(("Invalid block dimensions", kSize.width % 2 == 1 && kSize.height % 2 == 1 && kSize.width > 1 && kSize.height > 1));
        assert(("src and dst must be one channel grayscale images", src.channels() == 1, dst.channels() == 1));
        
        int width = src.cols;
        int height = src.rows;

        dst = cv::Mat1b( src.size() );
        
        cv::Mat1d means, sqrMeans, variances, avgVarianceMat;
        
        cv::boxFilter(src, means, CV_64F, kSize, cv::Point(-1,-1), true, cv::BORDER_REPLICATE);
        cv::sqrBoxFilter(src, sqrMeans, CV_64F, kSize, cv::Point(-1,-1), true, cv::BORDER_REPLICATE);

        cv::Mat1d means2 = means.mul(means);
        variances = sqrMeans - (means.mul(means));
        
        // estimate the noiseVariance
        cv::reduce(variances, avgVarianceMat, 1, cv::REDUCE_SUM, -1);
        cv::reduce(avgVarianceMat, avgVarianceMat, 0, cv::REDUCE_SUM, -1);
        noiseVariance = avgVarianceMat(0, 0) / (width * height);
        
        // perform wiener equation
        for (int r = 0; r < height; ++r){
            // get row pointers
            uchar const * const srcRow = src.ptr<uchar>(r);
            uchar * const dstRow = dst.ptr<uchar>(r);
            double * const varRow = variances.ptr<double>(r);
            double * const meanRow = means.ptr<double>(r);
            for (int c = 0; c < width; ++c) {
                dstRow[c] = cv::saturate_cast<uchar>(
                    meanRow[c] + cv::max(0., varRow[c] - noiseVariance) / cv::max(varRow[c], noiseVariance) * (srcRow[c] - meanRow[c])
                );
		    }
        }

    }  
    
}