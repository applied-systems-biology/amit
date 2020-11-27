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

#include "state_classification.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include "Segmentation.h"


namespace state_classification
{
    
    /**
     * Calculate the euklidian distance between two points.
     *
     * @param v point 1
     * @param w point 2
     * @param factor optional
     * @return euklidean distance of the two points
     */
    double euklideanDistance(cv::Point &v, cv::Point &w, std::optional<double> factor){

        /// distinguish for specified factor
        if ( factor.has_value() )
            return factor.value() * sqrt(pow(((double)(v.x - w.x)), 2.)+pow(((double)(v.y - w.y)), 2.));
        else
            return sqrt(pow(((double)(v.x - w.x)), 2.)+pow(((double)(v.y - w.y)), 2.));

    }

    /**
     *
     * @param region
     * @return average distance
     */
    double computeAverageDistanceFromGreenPixelsToCentroid(RegionP &region){

        double dist  = 0.;

        std::vector<cv::Vec2i> green = region.green_region_pixels;
        std::vector<cv::Vec2i>::iterator greenit;

        cv::Point c = region.getCentroidPoint();

        for(greenit = green.begin(); greenit != green.end(); greenit++){
            cv::Point g;
            g.y = greenit->val[0];
            g.x = greenit->val[1];
            double d = euklideanDistance(c, g);
            dist += d;
        }

        dist /= (double) green.size();

        return dist;

    }

    /**
     * Function computes several features for the all regions and saves them in the feature file. Additionally some of them are also saved in a sample
     * vector. But I don't see, that this one is used somewhere.
     * The features are mostly distances between centroids and/or pixels of the green, grey and overlay portions of the region as well as area and aspect ratio.
     *
     * 1. resize the two empty region vectors
     * 2. open feature output file
     * 3. go through all regions in regions
     *    - collect number of green and red pixels overlapping with the region
     *
     *    - if green pixels were found and the region is not too small
     * 		- clone current region and save number of green pixels
     * 		- compute centroid and save image number and coordinates of region centroid in feature file
     * 		- compute the centroid of the grey pixels of the region and save it in the feature file
     * 		- compute euklidean distance between the grey and overlay centroids and save it in the sample vector [0]
     * 		- compute centroid of only the green pixels and save it in the feature file
     * 		- compute euklidean distance between green and overlay centroids and save it in the sample vector [1]
     * 		- compute average distance from green pixels to overlay centroid and save it in the feature file and in the sample vector [2]
     * 		- save overlay area in feature file and sample vector [3]
     * 		- save overlay areas without green area in feature file and sample vector [4]
     * 		- save green area in feature file and sample vector [5]
     * 		- set ID of the region (increment ID)
     * 		- save ratio of green region pixels in feature file and sample vector [6]
     * 		- save ratio of red pixels in feature file
     * 		- fit ellipse around the region and save width and hight in feature file
     * 		- compute aspect ratio and save in feature file and sample vector [7]
     * 		- save ID in feature file
     * 		- set klass of the region to 2 and save the region in the fungi region vector
     *
     *    - if green == 0: save region in the PMN vector
     *
     *    --> all small regions will be missed with these two conditions and are therefore removed from the result
     *
     * 4. close feature output file
     * 
     * @param regions combined PMN and fungal regions
     * @param pmns empty region vector to save pmn regions
     * @param fungi empty region vector to save all regions overlapping with green regions
     * @param green fungal binary image
     * @param red red binary image
     * @param file output file of classification
     * @param t number of the current image +1
     * @param ID ID
     * @param f factor (depends on rows of the image) --> I think, this is used to normalize for different image sizes    
     */
    void computeFeatures(std::vector<Region> &regions, std::vector<Region> &pmns, std::vector<RegionP> &fungi, cv::Mat &green, cv::Mat &red, std::string &file, const int &t, int ID, const double &f){

        /// preprocessing
        pmns.resize(0);
        fungi.resize(0);

        std::ofstream outfile;
        outfile.open(file.c_str(), std::ofstream::app);

        std::vector<Region>::iterator regit;

        /// iterate over all regions of PMNs / immune cells and green fungi pathogens
        for(regit = regions.begin(); regit != regions.end(); regit++){
            std::vector<cv::Vec2i>::iterator pixit;

            int n_green = 0;
            int n_red = 0;
            /// iterate through all region pixels
            std::vector<cv::Vec2i> green_pixels;
            std::vector<cv::Vec2i> red_pixels;

            /// compute green and red pixels
            for(pixit = regit->region_pixels.begin(); pixit != regit->region_pixels.end(); pixit++){
                cv::Point p(pixit->val[1], pixit->val[0]);
                if(green.at<uchar>(p) > 0){
                    n_green++;
                    green_pixels.push_back(*pixit);
                }
                if(red.at<uchar>(p) > 0){
                    n_red++;
                    red_pixels.push_back(*pixit);
                }
            }


            if( n_green > 0 && regit->contour_pixels.size() > 5 && regit->region_pixels.size() > 10 ){

                std::vector<float> sample(8);
                RegionP newregion(*regit);

                /// set number of green pixels
                newregion.addGreenPxs(green_pixels);
                regit->computeCentroid();
                cv::Point v = regit->getCentroidPoint();

                outfile << t << "\t" << v.y << "\t" << v.x << "\t";

                newregion.computeCentroidFromGrayPixels();
                cv::Point w = newregion.getCentroidPoint();

                outfile << w.y << "\t" << w.x << "\t";
                outfile << euklideanDistance(v,w, f) << "\t";

                /// distance between overlay and brightfield centroids
                sample[0] = (float) euklideanDistance(v,w, f);

                newregion.computeCentroidFromGreenPixels();
                w = newregion.getCentroidPoint();

                outfile << w.y << "\t" << w.x << "\t";

                cv::Point u = newregion.getCentroidPoint();
                outfile << euklideanDistance(v, u, f) << "\t";

                /// distance between overlay and green centroids
                sample[1] = (float) euklideanDistance(v, u, f);

                /// average distance from green pixels to the centroid
                newregion.computeCentroid();
                outfile << computeAverageDistanceFromGreenPixelsToCentroid(newregion)*f << "\t";

                /// average distance of green pixels to centroid
                sample[2] = (float)  computeAverageDistanceFromGreenPixelsToCentroid(newregion)*f;

                ///// area /////

                outfile << newregion.getArea() * f * f<< "\t";
                sample[3] = newregion.getArea() * f * f; //complete area
                outfile <<   (newregion.getArea() - newregion.green_region_pixels.size()) * f * f << "\t";

                /// brightfield area
                sample[4] = (float) (newregion.getArea() - newregion.green_region_pixels.size()) * f * f;

                outfile <<  newregion.green_region_pixels.size() *f *f << "\t";

                /// green area
                sample[5] =  (float) newregion.green_region_pixels.size() *f *f;

                newregion.setId(++ID);
                outfile <<  newregion.getFungalRatioGreen() << "\t";

                /// ratio of green region pixels
                sample[6] = (float) newregion.getFungalRatioGreen() ;

                outfile <<  (double)  red_pixels.size() /(double)  newregion.getArea() << "\t";

                std::vector<cv::Point> points;
                for(size_t i = 0; i < regit->contour_pixels.size(); i++){
                    cv::Point p(regit->contour_pixels.at(i).val[1], regit->contour_pixels.at(i).val[0]);
                    points.push_back(p);
                }

                cv::RotatedRect r = cv::fitEllipse(points);
                outfile << r.size.height *f << "\t" << r.size.width *f << "\t";

                /// aspect ratio
                double ar = (double)r.size.width / (double) r.size.height ;
                outfile << ar << "\t";

                /// aspect ratio of enclosing ellipse
                sample[7] = (float) ar;

                /// shape factor
                double sf = (double) newregion.region_pixels.size() / pow((double)newregion.contour_pixels.size(), 2.);
                outfile << sf << "\t";

                outfile << newregion.getId() << std::endl;

                /// assign cluster cell
                newregion.setClass( cell_class::immune::CLUSTER );
                fungi.push_back(newregion);

            }
            /// remove noise
            else if(n_green == 0){
                pmns.push_back(*regit);
            }

        }

        outfile.close();

    }

    /**
     * Program to run the state classification part of the AMIT algorithm. 
     * (1) Read the file paths from the binary folder and from the green folder.
     * (2) Read images from the binary, fungal and dead cell directories.
     * (3) Perform state classification using binary, green and red segmented images.
     *
     * 1. create features output file and write header line
     * 2. read images from the three directories
     * 3. go through all images
     *    - add binary PMN images and fungal images into one image
     *    - do a morphological closing operation
     *    - collect Regions from the added image
     *    - compute features
     *    - collect number of green, pmn and other regions (but I don't see this numbers used somewhere?)
     *
     * 4. close feature output file
     * 5. generate classification output file
     * 6. perform classification with R using the feature file and save output in the classification file
     * 
     * @param images_binary images containing binary images of PMNs
     * @param images_fungi images containing binary images of fungal cells
     * @param images_dead images containing binary images of dead cells
     * @param START start number of images which will be taken into account
     * @param END end number of images which will be taken into account
     * @return path to 'features.csv'
     */
    std::string stateClassification(std::vector<cv::Mat> images_binary, std::vector<cv::Mat> images_fungi, std::vector<cv::Mat> images_dead, std::string &OUTPUTDIR,  const int &N_THREADS, const bool &FLAG_DEBUG){

        /**********************************************************/
        ///              calculate and store features              /
        /**********************************************************/

        std::stringstream ss;
        ss << OUTPUTDIR << "features.csv";

        std::string file = ss.str();
        std::ofstream outfile(file);

        /// write headline
        if(outfile.is_open()){
            std::stringstream ss_headline;
            ss_headline << "t" << "\t" << "cx" << "\t" << "cy" << "\t" << "cx_gray" << "\t" << "cy_gray" << "\t" << "d_gray" << "\t"
                << "cx_green" << "\t" << "cy_green" << "\t" << "d_green" << "\t" << "d_green2" << "\t" << "A" << "\t"
                << "A_gray" << "\t" << "A_green" << "\t" << "FR" << "\t" << "RR" << "\t" << "height" << "\t" << "width" << "\t"
                << "aspect_ratio" << "\t" << "sf" << "\t" << "ID";
            outfile << ss_headline.str() << std::endl;
        }
        outfile.close();

        const int ROWS = images_binary[0].rows;
        // const int COLS = images_binary[0].cols;

        cv::Mat samples;

        int n_green = 0, n = 0, n_pmn = 0;
        int ID = 0;

        // multi-threaded processing of brightfield images // TODO omp implementieren
        #pragma omp parallel for firstprivate(images_binary, images_fungi, n, n_green, n_pmn) num_threads(N_THREADS)
        for(size_t i = 0; i < images_binary.size(); i++){
            
            std::vector<Region> tmp, pmns;
            std::vector<RegionP> cgs, tmp_cgs;

            cv::Mat un;
            cv::add(images_binary[i], images_fungi[i], un);

            cv::morphologyEx(un, un, cv::MORPH_CLOSE, cv::Mat(), cv::Point(-1,-1),1);

            Segmentation::segmentRegionPixelBased2DRecursive(un, tmp);

            /// AMITv1: works only for quadratic images
            double FACTOR = 2048 / ROWS;

            state_classification::computeFeatures(tmp, pmns, cgs, images_fungi[i], images_dead[i], file, i+1, ID, FACTOR);

            n_green += cgs.size();
            n += cgs.size();
            n += pmns.size();
            n_pmn += pmns.size();

        }

        outfile.close();

        std::cout << " [2.2] Write features.csv to:\t" << OUTPUTDIR << std::endl;

        
        /// R - command
        std::string resultfile = OUTPUTDIR + "classification.csv";

        std::stringstream stringStream;
        stringStream << "Rscript ./src/classification.R " << file << " " << resultfile;
        
        std::string cmd = stringStream.str();

        std::cout << " [2.3] R-Command:\t" << cmd << std::endl;      

        try {

            int tmp = std::system(cmd.c_str());
            std::cout << " [2.4] std::system(cmd.c_str()):\t" << tmp << std::endl;
					     
        }
        catch(const std::exception& e) {
            std::cout << " [ERROR] by perform R-script" << std::endl;
            std::cerr << e.what() << std::endl;
        }       
        
        return file;

    }

}


