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

#include "segment_via_GMM.h"
#include <assert.h>
#include <stdlib.h>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include "IOputs.h"
#include "Segmentation.h"


namespace gmm
{

    /**
     * This function produces a binary image from a set of regions.
     * It creates a black image and adds the regions to the image with a gray value of 255.
     *
     * @param rows number of rows of the image
     * @param cols number of columns of the image
     * @param regions set of regions in the image
     *
     * @return binary image
     */
    cv::Mat create_binary_image_from_regions(const int &rows, const int &cols, std::vector<Region> &regions){
        cv::Mat im(rows, cols, CV_8UC1, cv::Scalar(0));

        std::vector<Region>::iterator reg_it;
        for(reg_it = regions.begin(); reg_it != regions.end(); reg_it++){
            std::vector<cv::Point> points;
            for(size_t ci = 0; ci < (size_t) reg_it->region_pixels.size(); ci++){
                cv::Point p(reg_it->region_pixels.at(ci)[1], reg_it->region_pixels.at(ci)[0]);
                im.at<uchar>(p) = 255;
            }
        }
        return im;
    }

    /**
     * Computes the temporal variance of the passed images.
     * First, it checks, if the number of images is uneven (which it has to be).
     * Then it extracts at each image position the pixel values of all images, 
     * computes mean and standard deviation and converts the std to variance. 
     * The variance is saved at the same image position in the result image.
     *
     * @param images vector (number of images used depends on the parameter N_TEMPVAR)
     * 
     * @return result image that contains the temporal variances at each pixel
     */
    cv::Mat compute_temporal_variance(const std::vector<cv::Mat> &images){
        int rows = images[0].rows;
        int cols = images[0].cols;

        int  t = images.size();
        // check for number of images even
        if(t % 2 == 0)
            std::cout << " [1.1.2] Error: imagestack for calculating temporal variance not of uneven size!" << std::endl;

        cv::Mat result = cv::Mat::zeros(rows, cols, CV_32FC1);
        
        #pragma omp parallel for shared(images)
        for(size_t i = 0; i < (size_t)rows; ++i){
            for(size_t j = 0; j < (size_t)cols; ++j){
                // essentially, a vector that contains pixel values of all images at position (i,j)
                cv::Mat tmp(t, 1, CV_8UC1, cv::Scalar(0)); 
                
                for(size_t tmp_i = 0; tmp_i < (size_t)t; ++tmp_i)
                    tmp.at<uchar>(tmp_i) = images[tmp_i].at<uchar>(i,j);
                
                cv::Mat mean, std;
                cv::meanStdDev(tmp, mean, std);
                
                // saturate_cast is used to round value to nearest value of the target data type
                // important to avoid introducing artefacts by rounding/casting errors into the image
                float value = cv::saturate_cast<float>(std.at<double>(0));
                value = pow(value, 2.);
                result.at<float>(i,j) = value;
            }
        }

        return result;
    }

    /**
     * This function computes the local/spatial image variance.
     * The variance is calculated within a window of size kerneldim x kerneldim.
     * Image will convert to float32, calculate, mean, sigma, variance.
     *
     * @param image input image
     * @param kerneldim dimension of the kernel
     * 
     * @return variance the spatial variance image of the input
     */
    cv::Mat compute_local_variance(const cv::Mat &image, const int &kerneldim){	
        cv::Mat image32f;
        image.convertTo(image32f, CV_32F);

        cv::Mat mu;
        cv::blur(image32f, mu, cv::Size(kerneldim, kerneldim));
        cv::Mat mu2;
        cv::blur(image32f.mul(image32f), mu2, cv::Size(kerneldim, kerneldim));
        cv::Mat sigma;
        cv::sqrt(mu2 - mu.mul(mu), sigma);
        
        cv::Mat variance;
        cv::pow(sigma, 2., variance);
        
        return variance;
    }

    /**
     * Classifies pixels in the current image according to their temporal and spatial variance and
     * uses its own em_model and label objects
     * 
     * @param n_i image index
     * @param res result image
     * 
     */
    void classify_pixel(const int &n_i, cv::Mat &dst, const std::vector<cv::Mat> &images, cv::Ptr<cv::ml::EM> &em_model, const int &N_TEMPVAR, const int &C, const bool &FLAG_DEBUG) {
        std::vector<cv::Mat> imagestack;
        cv::Mat image = images[n_i].clone();
        
        assert (!image.empty());
        
        // set "neighboring" images in array 
        for(int index = -(N_TEMPVAR/2); index <= (N_TEMPVAR/2); ++index){
            imagestack.push_back(images[n_i + index]);
            if (images[n_i+index].rows <= 0 or images[n_i+index].cols <= 0)
                std::cout << " [1.1.1] Colour image doesn't exist: image " << n_i << ", temporal variance image " << index << std::endl;
        }
        
        // compute temporal variance      
        cv::Mat imtempvar = gmm::compute_temporal_variance(imagestack);    
        
        // compute spatial/local variance
        cv::Mat imvar = gmm::compute_local_variance(image, 3); 
        
        if(FLAG_DEBUG)
            io::draw_images(imtempvar, imvar);
        
        // compute samples for learning the GMM (flatten)
        int n_samples = (image.cols-2)*(image.rows-2);
        cv::Mat samples = cv::Mat::zeros( n_samples, 2, CV_32FC1 );
        
        int k = 0;
        for(int i = 1; i < image.rows-1; ++i){
            for(int j = 1; j < image.cols-1; ++j){                
                samples.at<float>(k,0) = imvar.at<float>(i,j);
                samples.at<float>(k,1) = imtempvar.at<float>(i,j);
                
                ++k;            
            }
        }
                
        // declare GMM (cv::ml::EM) for temp/spatial variation
        em_model = cv::ml::EM::create();
        em_model->setClustersNumber(C); 
        em_model->setCovarianceMatrixType(cv::ml::EM::COV_MAT_GENERIC);
        em_model->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER+cv::TermCriteria::EPS, 300, 0.1));
            
        // output of EM-model
        cv::Mat labels;
        cv::Mat probs;

        // learn gmm-model
        em_model->trainEM( samples, cv::noArray(), labels, probs );

        //output of gmm parameters
		std::cout << " [1.2] GMM-means: " << em_model->getMeans() << std::endl;
		std::cout << " [1.3] GMM-weights: " << em_model->getWeights() << std::endl;

        cv::Mat result = cv::Mat::zeros(image.rows, image.cols, CV_8UC1);
                
        // determine if labels are in the right order (0,1,2,...)
        cv::Mat weights = em_model->getWeights();
        cv::Mat weights_sort;
        cv::sort(weights, weights_sort, cv::SortFlags::SORT_EVERY_ROW | cv::SortFlags::SORT_DESCENDING);
                
        // check if classes are in an ascending order
        cv::Mat weights_sub;
        cv::subtract(weights, weights_sort, weights_sub);
        weights_sub = cv::abs(weights_sub);
        
        double sum = 0;
        for(size_t i = 0; i < (size_t)C; ++i)
            sum += weights_sub.at<double>(0,i);

        // change labeling if it is not in an increasing order
        if(sum > 0.001){
            std::vector<int> indizes(C,0); //lookuptable of indizes

            //fill lookup-table
            for(size_t i = 0; i < (size_t)C; ++i){ 
                double temp = weights.at<double>(0,i);
                for(size_t j = 0; j < (size_t)C; ++j){ 
                    if(abs(temp - weights_sort.at<double>(j)) < 0.001){
                        indizes.at(i) = j;
                    }
                }
            }
            em_model->getWeights() = weights_sort;

            // change rows in means matrix
            cv::Mat means = em_model->getMeans();
            cv::Mat means_sort;
            means.copyTo(means_sort);

            for(size_t i = 0; i < (size_t)C; ++i){ 
                for(size_t j = 0; j < (size_t)C; ++j){
                    means_sort.at<double>(i, j) = means.at<double>(indizes.at(i),j);
                }
            }

            em_model->getMeans() = means_sort;

            for(size_t i = 0; i < (size_t)labels.rows; ++i){
                int j = labels.at<int>(i);
                labels.at<int>(i) = indizes.at(j);
            }
        }
        
        // output generation --> different colors for the different classes
        k = 0;
        for(size_t i = 1; i < (size_t)image.rows-1; ++i){
            for(size_t j = 1; j < (size_t)image.cols-1; ++j){
                int l = labels.at<int>(++k);
                if(l == 0){
                    result.at<uchar>(i,j) = 0;
                }
                else{
                    uchar v = 256/(C-1)*l; 
                    result.at<uchar>(i,j) = v-1;
                }


            }
        }
        
        result.copyTo(dst);

    }

    /**
     * Calls a segmentation function from the Segmentation library to segment regions in the image.
     * Afterwards all regions are checked. Only regions larger than 20 and not line-shaped are kept and saved in the passed region vector.
     *
     * @param image input image
     * @param regions vector of regions
     */
    void segment_PMNs_pixel_based_2D(const cv::Mat &image, std::vector<Region> &regions){

        std::vector<Region> tmp;
        Segmentation::segmentRegionPixelBased2DRecursive(image, tmp);

        for(std::vector<Region>::iterator regit = tmp.begin(); regit != tmp.end(); regit++){
            if(regit->getArea() > 20 && !gmm::is_line(*regit)){
                regions.push_back(*regit);
            }
        }

    }

    /**
     * This function determines wether a region is line-shaped or not.
     * Therefore, an enclosing rotated rectangle is computed that incorporates the region.
     * The aspect ratio calculated by longer side / shorter side (which is always < 1).
     * If the aspect ration is smaller .25 or the shorter side is less than 10 px, the function returns true.
     *
     * @param region region
     * @return true, if the region is line-shaped
     */
    bool is_line(const Region &region){
        std::vector<cv::Point> points;

        for(size_t ci = 0; ci < (size_t)region.region_pixels.size(); ci++){
            cv::Point p(region.region_pixels.at(ci)[1], region.region_pixels.at(ci)[0]);
            points.push_back(p);
        }


        // try-catch-clause added by Philipp Praetorius to prevent
        //  OpenCV Error: Incorrect size of input array (There should be at least 5 points to fit the ellipse) in fitEllipse
        cv::RotatedRect rect;
        try
        {
            rect = cv::fitEllipse(points);
        }
        catch(const std::exception& e)
        {
            std::cout << "cv::fitEllipse: There should be at least 5 points to fit the ellipse in function fitEllipse: " << points.size() << " Points" << std::endl;
            std::cerr << e.what() << '\n';

            std::cout << "CAUTION: gmm::is_line return false" << std::endl;
            return false;                
        }


        float ar;
        float dmin;

        if(rect.size.height < rect.size.width ){
            ar = rect.size.height / rect.size.width;
            dmin = rect.size.height;
        }
        else{
            ar = rect.size.width / rect.size.height;
            dmin = rect.size.width;
        }

        if(ar < .25){
            return true;
        }
        else if(ar < .75 && dmin < 10){
            return true;
        }
        return false;
    }

    /**
     * Function takes an image from the classification step and extracts PMNs by thresholding for the third class.
     * The PMN regions are closed with a kernel of size 3x3 for N_CLOSINGS number of times (probably to smooth out bumps in the regions?).
     * The resulting objects are extracted to a vector of Regions. In these Regions, remaining holes are filled and a binary image containing them is created.
     * The result is saved in the second input parameter with gray values of 255 for the PMN pixels.
     *
     * @param image with 3 different pixel values presenting the three classes
     * @param result binary image containing the PMN regions
     *
     */
    void detect_PMNs(const cv::Mat &src, cv::Mat &dst, const int &N_CLOSINGS, const bool &clear_border, const bool &FLAG_DEBUG){
        cv::Mat img_thresh, img_close;
        
        // threshold image to get only px that are classified to class 3
        cv::threshold(src, img_thresh, 200, 255, cv::THRESH_BINARY);
        
        // apply morphological closing (original version from Susanne)
        cv::Mat kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(3,3));
        
        cv::morphologyEx(img_thresh, img_close, cv::MORPH_CLOSE, kernel, cv::Point(-1,-1), N_CLOSINGS);
        
        // Test for other closing operation using variable kernel and not several operations ...
        
        //visual::createGUI(im, false, visual::imgHandling::none);

        if(FLAG_DEBUG){
            cv::imshow("img_thresh", img_thresh);
            cv::imshow("img_close", img_close);
            cv::waitKey(0);
        }

        std::vector<Region> regions;
        // segment regions from binary images
        gmm::segment_PMNs_pixel_based_2D(img_close, regions);

        // fill holes inside regions
        Segmentation::fill_holes(regions, clear_border);

        // create binary images from filled regions
        dst = create_binary_image_from_regions(img_close.rows, img_close.cols, regions);
    } 
 
    /**
     * Main-Function for GMM-segmentation
     *
     * @param tid thread ID
     * @param inDir input directory path
     * @param outDir output directory path
     *
     * - reads image at the index of the thread id
     * - calculates which images have to be processed by the current thread
     * - calculates which images are excluded because of the temporal variance
     * - goes through the calculated images
     * - do gmm classification on the image and save result in gmm image
     * - detect PMNs and save the in the binary image
     * - save the binary image
     */
    void gmm_segmentation(const std::vector<cv::Mat> &images, const int &N_THREADS, const int &N_TEMPVAR, const int &C, const int &N_CLOSINGS, const std::string &OUT_BINARY, const bool &clear_border, const bool &FLAG_DEBUG) {
                
        int first = N_TEMPVAR/2;
        int last = (signed) images.size() - first - 1;

        // multi-threaded processing of brightfield images
        #pragma omp parallel for firstprivate(images) num_threads(N_THREADS)
        for (size_t i = 0; i < images.size(); ++i)
        {
            // compute gmm of the second to penultimate image
            if(i >= (size_t)first && i <= (size_t)last) {
                cv::Mat img_gmm;

                // classify image pixels into three classes
                cv::Ptr<cv::ml::EM> em_model;
                gmm::classify_pixel(i, img_gmm, images, em_model, N_TEMPVAR, C, FLAG_DEBUG); 

                if(FLAG_DEBUG){
                    cv::imshow("img_gmm", img_gmm);
                    cv::waitKey(0);
                }
                                
                // detect PMNs in three-class-image
                cv::Mat binary;
                gmm::detect_PMNs(img_gmm, binary, N_CLOSINGS, clear_border, FLAG_DEBUG);
                                
                // save binary image including PMNs
                io::write_image(OUT_BINARY, binary, i, true);
            }
        }            
        
    }      
    
    
} // gmm
