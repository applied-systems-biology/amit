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

#pragma once

#include <vector>
#include <opencv2/opencv.hpp>


namespace cell_class 
{
    
    enum class immune { 
        NOISE = 0,
        SINGLE = 1, 
        CLUSTER = 2 
    };

    enum class interaction {
        FUNGAL_CELLS = 1,
        IMMUNE_TOUCHING_PATHOGEN = 2,
        PHAGOCYTOSED = 3,
        PHAGOCYTOSED_INTERACTION_IMMUNE = 4,
        PHAGOCYTOSED_INTERACTION_PATHOGEN = 5,
        DEAD_CELLS = 6
    };

    std::string getCellClass(immune c);
    immune setCellType(int c);
    std::string getCellClassP(interaction c);
    interaction setCellTypeP(int c);

};


class Region
{
    
    protected:

        int ID;                             /**< ID */
        cell_class::immune klass;           /**< class: 0 -> single cell, 1 -> cell cluster, 2 -> other */
        cv::Vec2f centroid;                 /**< centroid */            
                
        // bool interaction;                /**< 0 -> no, 1 -> yes */
        std::vector<int> interactionIDs;    /**< vector that saves the IDs of the cell track(s) with that the cell is currently interacting */
        cv::Mat im;                         /**< binary minimal image of the region */
        bool flat = false;                  /**< flag indicates if the region represent a flat cell */
        float minAxe = NAN;                 /**< ellipsoid major axis of the region */
	    float majAxe = NAN;                 /**< ellipsoid minor axis of the region */
        
    public:
        // constructor
        Region(/* args */)
        {
            this->ID = 0;
            /// single cells is the default class for single cells
            this->klass = cell_class::immune::SINGLE;
        }

        // destructor
        ~Region() {};

        // public: class - attributes
        int n_regions;
        std::vector<double> likelihood;
        std::vector<cv::Vec2i> region_pixels; /**< set of pixels defining the region */
        std::vector<cv::Vec2i> contour_pixels; /**< set of pixels defining the contour*/

        cv::Scalar color_border; /**< color of region border */

        // public: properties / class - methods
        cell_class::immune getClass();
        int getId();
        int getArea();
        cv::Mat getImage(); 
        cv::Scalar getCentroidScalar();
        cv::Point getCentroidPoint();
        void getMinMax(int &minx, int &miny, int &maxx, int &maxy);   
        void getPrincipalAxesMinEnclosingEllipse(cv::Size2f &s);   
        bool isInteracting();
        void getInteractionIds(std::vector<int> &ids);
        float getMaj();
	    float getMin();
	    bool getFlat();
        
        void addRegionPixels(const std::vector<cv::Vec2i> &v);
        void addContourPixels(const std::vector<cv::Vec2i> &v);
        void addRegionPixel(const cv::Vec2i &p);
        void addInteractionId(const int &id);

        bool hasNeighbor(cv::Vec2i &v, std::vector<cv::Vec2i> &list);
        friend int overlapN(Region &r, Region &p);
        friend int overlap(Region &r, Region &p);
        friend int overlap(Region *r, Region *p);

        void setClass(cell_class::immune c);
        void setId(const int &id);
        void setCentroid(const double &x, const double &y);
        void setFlat(const bool f);
        void setMinMajAxe(const float min, const float maj);
        
        void computeContour();
        void computeCentroid();
        friend double computeDistance(Region &p, Region &q);
        friend double computeDistance(Region *p, Region *q);
        void createImage();
        void changeInteractionID(int id, int id_new);        

};
