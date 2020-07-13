/*  
Copyright by Jan-Philipp_Praetorius

Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo
Figge
https://www.leibniz-hki.de/en/applied-systems-biology.html
HKI-Center for Systems Biology of Infection
Leibniz Institute for Natural Product Research and Infection Biology -
Hans Knöll Insitute (HKI)
Adolf-Reichwein-Straße 23, 07745 Jena, Germany 
 */

#include "Segmentation.h"


namespace Segmentation
{

    /**
     * Clockwise region growing approach (iterative approach).
     * This function is more expensive in time than the recursive approach but there is no limited recursion depth.
     *
     * @param im image
     * @param i current x-position
     * @param j current y-position
     * @param pixels set of foreground region pixels
     */
    void getRegion2DIterative(cv::Mat &im, int i, int j, std::vector<cv::Vec2i> &pixels){
        //add current pixel to the region
        cv::Vec2i v(i,j);
        pixels.push_back(v);

        for(size_t k = 0; k < (size_t) pixels.size(); k++){
            i = pixels.at(k).val[0];
            j = pixels.at(k).val[1];
                v.val[0] = i;
                v.val[1] = j;

                im.at<uchar>(i,j) = 200;

                //right
                if(j < im.cols-1 && im.at<uchar>(i,j+1) == 255){
                    im.at<uchar>(i,j+1) = 127;
                    cv::Vec2i v1(i,j+1);
                    pixels.push_back(v1);
                }
                //below
                if((i < im.rows - 1) && im.at<uchar>(i+1, j) == 255){
                    im.at<uchar>(i+1,j) = 127;
                    cv::Vec2i v1(i+1,j);
                    pixels.push_back(v1);
                }
                //left
                if((j > 0) && im.at<uchar>(i,j-1) == 255){
                    im.at<uchar>(i,j-1) = 127;
                    cv::Vec2i v1(i,j-1);
                    pixels.push_back(v1);
                }
                //above
                if((i > 0) && im.at<uchar>(i-1, j) == 255){
                    im.at<uchar>(i-1,j) = 127;
                    cv::Vec2i v1(i-1,j);
                    pixels.push_back(v1);
                }
        }
    }

     /**
     * Clockwise region growing approach.
     *
     * @param im image
     * @param i current x-position
     * @param j current y-position
     * @param pixels set of foreground region pixels
     */
    void getRegion2DRecursive(cv::Mat &im, int i, int j, std::vector<cv::Vec2i> &pixels){
        // set current pixel to gray
        im.at<uchar>(i,j) = 127;
    
        // add current pixel to the region
        cv::Vec2i v(i,j);
        pixels.push_back(v);

        // right
        if((j < im.cols - 1) && im.at<uchar>(i,j+1) == 255){
            getRegion2DRecursive(im, i, j+1, pixels);
        }
        // below
        if((i < im.rows - 1) && im.at<uchar>(i+1, j) == 255){
            getRegion2DRecursive(im, i+1, j, pixels);
        }
        // left
        if((j > 0) && im.at<uchar>(i,j-1) == 255){
            getRegion2DRecursive(im, i, j-1, pixels);
        }
        // above
        if((i > 0) && im.at<uchar>(i-1, j) == 255){
            getRegion2DRecursive(im, i-1, j, pixels);
        }
    }

    /**
     * This function estimates the contour pixels of a region given a binary image.
     * Contour pixels have at least one direct neighbor that is black.
     * time complexity: O(#region pixels)
     *
     * @param image binary image
     * @param region_pixels pixels of current region
     * @param pixels contour pixels
     */
    void getContour2D(cv::Mat &image, const std::vector<cv::Vec2i> &region_pixels, std::vector<cv::Vec2i> &pixels){
        cv::Mat im;
        image.copyTo(im);

        std::vector<cv::Vec2i>::const_iterator it;
        for(it = region_pixels.begin(); it != region_pixels.end(); it++){
            int i = (*it)[0];
            int j = (*it)[1];

            // pixel at image border?
            if((j == im.cols - 1) || (i == im.rows -1) || (i == 0) || (j == 0) ){
                pixels.push_back(*it);
                im.at<uchar>(i,j) = 77;
            }
            else{// 8-neighborhood (Moore)
                cv::Vec2i v(i,j);
                if((im.at<uchar>(i,j+1) == 0)){
                    pixels.push_back(v);
                    im.at<uchar>(i,j) = 77;
                }
                else if(im.at<uchar>(i+1,j+1) == 0){
                    pixels.push_back(v);
                    im.at<uchar>(i,j) = 77;
                }
                else if(im.at<uchar>(i+1,j-1) == 0){
                    pixels.push_back(v);
                    im.at<uchar>(i,j) = 77;
                }
                else if(im.at<uchar>(i+1,j) == 0){
                    pixels.push_back(v);
                    im.at<uchar>(i,j) = 77;
                }
                else if(im.at<uchar>(i,j-1) == 0){
                    pixels.push_back(v);
                    im.at<uchar>(i,j) = 77;
                }
                else if(im.at<uchar>(i-1,j+1) == 0){
                    pixels.push_back(v);
                    im.at<uchar>(i,j) = 77;
                }
                else if(im.at<uchar>(i-1,j) == 0){
                    pixels.push_back(v);
                    im.at<uchar>(i,j) = 77;
                }
                else if(im.at<uchar>(i-1,j-1) == 0){
                    pixels.push_back(v);
                    im.at<uchar>(i,j) =77;
                }
            }
        }
    }

    /**
     * returns true or false whether the vector of pixels contains a pixel at the image borders
     */
    bool touchesBorder(const std::vector<cv::Vec2i> &pixels, const cv::Mat &image){
        int rows = image.rows-1;
        int cols = image.cols-1;

        std::vector<cv::Vec2i>::const_iterator it;
        for(it = pixels.begin(); it != pixels.end(); it++){
            if(it->val[0] == 0 || it->val[1] == 0 || it->val[0] == rows || it->val[1] == cols){
                return true;
            }
        }
        return false;
    }    

    /**
     * This function fills holes within a set of regions
     *
     * @param regions set of regions
     */
    void fill_holes(std::vector<Region> & regions, const bool &clear_border){
        std::vector<Region>::iterator reg_it;

        for(reg_it = regions.begin(); reg_it != regions.end(); reg_it++){
            reg_it->createImage();
            cv::Mat im = reg_it->getImage();
            // compute inverse image --> foreground pixels become background and vice versa
            cv::bitwise_not(im,im); 

            // now holes are detected as foreground objects
            std::vector<Region> holes;
            Segmentation::segmentRegionPixelBased2DIterative(im, holes, clear_border);

            if(holes.size() > 0){
                int minx, miny, maxx, maxy;
                reg_it->getMinMax(minx, miny, maxx, maxy);
    
                // add hole pixels to the region
                std::vector<Region>::iterator reg_it2;
                for(reg_it2 = holes.begin(); reg_it2 != holes.end(); reg_it2++){
                    for(size_t ci = 0; ci < (size_t) reg_it2->region_pixels.size(); ci++){
                        cv::Vec2i v(reg_it2->region_pixels.at(ci)[0] + minx, reg_it2->region_pixels.at(ci)[1]+ miny);
                        reg_it->addRegionPixel(v);
                    }
                }

                // recompute contour
                reg_it->computeContour(); 
                reg_it->createImage();
            }
        }
    }

    /**
     * This function performs the segmentation of a binary image  using a region-growing approach.
     * The set of regions is returned by a call-by-reference of the parameter regions.
     *
     * @param image image
     * @param regions set of regions (call-by-reference)
     * @param clear_border removes objects which are connected to the image border
     */
    void segmentRegionPixelBased2DIterative(const cv::Mat &image, std::vector<Region> &regions, const bool &clear_border){
        // initialize
        regions.resize(0);

        // hard copy
        cv::Mat im = image.clone(); 
        cv::Mat im2 = image.clone();

        int rows = im.rows;
        int cols = im.cols;

        // get regions and contours and save them as Region2D objects
        if(im.isContinuous()){
            cols *= rows;
            rows = 1;
        }

        for(size_t i = 0; i < (size_t) im.rows; i++){
            uchar *data = im.ptr(i);

            for(int j = 0; j < (size_t) im.cols; j++){
                std::vector<cv::Vec2i> pixels_c, pixels_r;

                uchar value = *data;
                // search for white pixels
                if(value == 255){

                    Region reg;
                    
                    Segmentation::getRegion2DIterative(im, i, j, pixels_r);
                    Segmentation::getContour2D(im2, pixels_r, pixels_c);

                    /// discard objects at the image border if choosen
                    if (clear_border) {
                        if (Segmentation::touchesBorder(pixels_c, im2)) {
                            continue;
                        }
                    }

                    reg.addRegionPixels(pixels_r);
                    reg.addContourPixels(pixels_c);
    
                    pixels_r.resize(0);
                    pixels_c.resize(0);

                    reg.computeCentroid();
                    reg.createImage();


                    regions.push_back(reg);    
                }

                data++;
            }
        }

    }

    /** // TODO evtl ganze Funktion löschen ... nicht verwendet
     * This function performs the segmentation of a binary image  using a region-growing approach.
     * The set of regions is returned by a call-by-reference of the parameter regions.
     *
     * @param image image
     * @param regions set of regions (call-by-reference)
     */
    void segmentRegionPixelBased2DRecursive(const cv::Mat &image, std::vector<Region> &regions){
        // initialize
        regions.resize(0);

        // hard copy
        cv::Mat im = image.clone(); 
        cv::Mat im2 = image.clone();

        int rows = im.rows;
        int cols = im.cols;

        // get regions and contours and save them as Region2D objects
        if(im.isContinuous()){
            cols *= rows;
            rows = 1;
        }

        for(size_t i = 0; i < (size_t)im.rows; i++){
            uchar *data = im.ptr(i);

            for(size_t j = 0; j < (size_t)im.cols; j++){
                std::vector<cv::Vec2i> pixels_c, pixels_r;

                uchar value = *data;
                // search for white pixels
                if(value == 255){
                    Region reg;
                    Segmentation::getRegion2DRecursive(im, i, j, pixels_r);
                    Segmentation::getContour2D(im2, pixels_r, pixels_c);

                    if(Segmentation::touchesBorder(pixels_c, im2))
                        continue;

                    reg.addRegionPixels(pixels_r);
                    reg.addContourPixels(pixels_c);
    //				reg.computeContour(); 
                    // TODO: Philipp: alles was hier aussortiert ist auf 
                    // Notwendigkeit prüfen und ggf. löschen

                    pixels_r.resize(0);
                    pixels_c.resize(0);

                    reg.computeCentroid();
                    reg.createImage();
    //				reg.showImage();

    //				for(unsigned int i = 0; i < reg.region_pixels.size(); i++){
    //					cout << reg.region_pixels.at(i)[0] << " " << reg.region_pixels.at(i)[1] << endl;
    //				}
                    //sort pixels (erst nach x und dann nach y)
    //				sort(reg.region_pixels.begin(), reg.region_pixels.end(), Segmentation::sortVec2i);
    //				for(unsigned int i = 0; i < reg.region_pixels.size(); i++){
    //					cout << reg.region_pixels.at(i)[0] << " " << reg.region_pixels.at(i)[1] << endl;
    //				}

    //				if(!Segmentation::touchesBorder(reg.contour_pixels, im)){
                    regions.push_back(reg);
    //				}

                }
                data++;
            }
        }

    }
     
} // Segmentation
 