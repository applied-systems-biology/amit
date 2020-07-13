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

#pragma once

#include <vector>
#include <algorithm>  
#include <opencv2/opencv.hpp>
#include <numeric>
#include "Region.h"


/**
 * Contains the information for one image. 
 * Each element one of the vectors represents one region/area within the image.
 */
struct feature_data
{
    std::vector<std::vector<cv::Point>> PixelIdxList;  // list of pixels for each object on binary image
    std::vector<int> ID;                               // current object ID
    std::vector<std::vector<int>> IDin;                // vector of object(s) ID(s) which compose current object
    std::vector<int> IDout;                            // vector of object(s) ID(s) which came out of current object

    std::vector<int> In;                               // number of objects, which comes in current object
    std::vector<int> Out;                              // number of objects, which comes out of current object
    std::vector<int> t;                                // current time point
    std::vector<int> type;                             // type of object (1: single object, 2: cluster)
    std::vector<int> num_regions;                      // number of regions within a single ROI/object

    /// constructors
    feature_data() {}
    feature_data(const std::vector<std::vector<cv::Point>> &pxl_idx_list) : PixelIdxList(pxl_idx_list) {} 
   
    /**
     * Append all of a given object with all detailed lists to the existing one.
     *      
     * @param objList list which should be appended to the existing one.
     */
    void append(const feature_data &objList){
        this->PixelIdxList.insert( this->PixelIdxList.end(), objList.PixelIdxList.begin(), objList.PixelIdxList.end() );
            
        this->ID.insert( this->ID.end(), objList.ID.begin(), objList.ID.end() );
        this->IDin.insert( this->IDin.end(), objList.IDin.begin(), objList.IDin.end() );
        this->IDout.insert( this->IDout.end(), objList.IDout.begin(), objList.IDout.end() );

        this->In.insert( this->In.end(), objList.In.begin(), objList.In.end() );
        this->Out.insert( this->Out.end(), objList.Out.begin(), objList.Out.end() );
        this->t.insert( this->t.end(), objList.t.begin(), objList.t.end() );
        this->num_regions.insert( this->num_regions.end(), objList.num_regions.begin(), objList.num_regions.end() );
    }

    /**
     * Determine the number of regions per cluster. 
     * 
     * For special cases, e.g. cell leave one cluster to merge directly with another cluster, apply a heuristic. 
     * 
     * (1) collect the size of all single cells from forward projection 
     * (2) iterate over backwards-projection to detect size of splitting objects in forward-projection
     * (3) iterate over all objects in objList-fw and bw by time to get all clusters already exist in the first frame
     * (4) iterate all entries and enter the number of regions within a ROI/objects
     *      
     * @param objList_bw object list of backward-projection.
     */
    void getNumberRegions(feature_data &objList_bw){

        /// map with key = ID, value = number of regions
        std::map<int, int> id_numRegions;

        int t_max = *std::max_element(this->t.begin(), this->t.end() );

        // (1) collect the size of all single cells from forward projection 
        std::vector<double> cell_size;

        for (size_t i = 0; i < this->ID.size(); ++i) {
            /// check if cell is a single cells, if yes: collect size, skip to small objects (when they are at the image boundary)
            if(this->type[i] == 1){
                double area = fabs(contourArea(cv::Mat( objList_bw.PixelIdxList[i] )));

                if (area > 30.)
                    cell_size.push_back( area );
            }
        }

        double cell_average_size = std::accumulate( cell_size.begin(), cell_size.end(), 0.0)/cell_size.size();

        // (2) iterate over backwards-projection to detect size of splitting objects in forward-projection
        std::map<int, int> id_numRegions_bw;
        for (size_t j = 0; j < objList_bw.ID.size() ; ++j) {

            int id_bw = objList_bw.ID[j];

            /// ID-bw already exists
            if ( id_numRegions_bw.count( id_bw ) ){
                continue;
            }
            /// more than 1 incoming objects => generation of cluster
            else if (objList_bw.In[j] > 1) {
                int n_regions = 0;

                /// get for each incoming objects the corresponding number of regions (1 for single cell)
                for (int id_in : objList_bw.IDin[j]) {
                    n_regions += id_numRegions_bw[id_in];
                }

                id_numRegions_bw[ id_bw ] = n_regions;

            }
            /// cell just remains or more than 1 outgoing cell => split cluster
            else {
                id_numRegions_bw[objList_bw.IDout[j]] = objList_bw.num_regions[j];
            }
        }

        ////////// detect all cluster generations and splittings //////////
        // (3) iterate over all objects in objList-fw and bw by time to get all clusters already exist in the first frame
        for (int i = 0; i < (int)this->ID.size(); ++i) {

            int id = this->ID[i];

            int t_fw = this->t[i];
            int t_bw = t_max - t_fw;

            /// get index number where t_fw = t_bw in backwards projection (ATTENTION: strong dependency on index sequence)
            std::vector<int>::iterator it_bw = std::find(objList_bw.t.begin(), objList_bw.t.end(), t_bw);
            int idx_t_bw_first = std::distance( objList_bw.t.begin(), it_bw);

            /// calculate the same index in backwards projection as i is for forward-projection
            std::vector<int>::iterator it_fw = std::find(this->t.begin(), this->t.end(), t_fw);
            int idx_t_fw_first = std::distance( this->t.begin(), it_fw);

            int idx_t = abs( i - idx_t_fw_first );
            int i_bw = idx_t_bw_first + idx_t;


            /// 1. check if a cluster split occur (Out > 1) and a single cell remains (IDout: type = 1)
            if ( id != this->IDout[i] and this->Out[i] > 1 and this->type[i] == 1 ){
               id_numRegions[ this->IDout[i] ] = 1;
               continue;
            }
            /// 2. check if a cluster split occur (Out > 1) and a cluster remains (IDout: type = 2)
            else if (this->Out[i] > 1 and this->type[i] == 2 ){
                std::cout << "INFO: a cluster split occur and a cluster remains (IDout: type = 2) i: " << i << std::endl;

                /// CAUTION: special case when object jumps from one cluster to another

                /// assign generation cluster of backward projection to split-cluster in forward projection
                int num_regions_bw = id_numRegions_bw[ objList_bw.ID[i_bw] ];

                /// compute average size of remaining cluster region, round up at .5 to next integer
                double area = fabs(contourArea(cv::Mat( this->PixelIdxList[i] )));
                int num_regions_average = std::round( area / cell_average_size );

                /// ensure to have number inside cluster at least 2: because outgoing object is still a cluster
                std::vector<int> tmp_num_regions;
                if (num_regions_bw > 1) tmp_num_regions.push_back( num_regions_bw );
                if (num_regions_average > 1) tmp_num_regions.push_back( num_regions_average );

                if ( tmp_num_regions.size() == 0 ){
                    tmp_num_regions.push_back( 2 );
                }

                /// assign the min-number regions for the remaining cluster, when a region of one cluster jumps to another the
                /// error of incorrectly assigning numbers for individual cells is multiplied by the error for all later frames
                int num_min_regions = *std::min_element(tmp_num_regions.begin(), tmp_num_regions.end() );

                id_numRegions[ this->IDout[i] ] = num_min_regions;

                continue;
            }
            /// 3. check if incoming object(s) already is a cluster and add the number of regions accordingly
            if ( this->In[i] > 1 ){

               int n_regions = 0;

               /// get for each incoming objects the corresponding number of regions (1 for single cell)
               for( int id_in : this->IDin[i] ) {
                   n_regions += id_numRegions[ id_in ];
               }

               id_numRegions[ id ] = n_regions;

            }
            /// 4. check if ID is already stored, to avoid overwrite ID - info
            if (id_numRegions.count( id )) {
               continue;
            }

            /// 5. check if a cluster already exist in the first frame, by an overlay of forward and backwards list
            /// assign all of the other unique ID's there corresponding number of regions
            int num_regions_fw = this->In[ i ];
            int num_regions_bw = objList_bw.In[ i_bw ];
            std::vector<int> num_regions_tmp { num_regions_fw, num_regions_bw} ;

            int num_regions_max = *std::max_element(num_regions_tmp.begin(), num_regions_tmp.end() );

            /// store id with number of regions, when cluster is generated
            id_numRegions[ id ] = num_regions_max;

        }


        ////////// assign number of regions to certain ID's //////////
        // (4) iterate all entries and enter the number of regions within a ROI/objects (single cell == 1, otherwise > 1)
        for (size_t i = 0; i < this->ID.size(); ++i) {

            int id = this->ID[i];

            /// check for cluster and enter the number of regions
            if (id_numRegions.count( id )){

                /// check if an objects is splitted ( ID != IDout )
                if( id != this->IDout[i] ){
                    this->num_regions[i] = id_numRegions[ this->IDout[i] ];
                } else {
                    this->num_regions[i] = id_numRegions[id];
                }
            }
        }

   }

    /**
     * Get the index of cluster-regions and the number of regions within a cluster from a given time frame.
     *      
     * @param t specified time point / frame.
     * @param indices all detected indices of cluster per frame t.
     * @param num_regions the number of regions per time point index. 
     */
    void getClusterIndex(const int t, std::vector<int> &indices, std::vector<int> &num_regions){

        /// iterate over all regions in time = t
        for (size_t i = 0; i < this->t.size(); ++i) {
            /// check if frame is the given frame and the region is a cluster (=2)
            if (this->t[i] == t and this->type[i] == 2){
                indices.push_back(i);

                num_regions.push_back( this->num_regions[i] );
            }
        }

    }

};


namespace ibp_cluster_detection
{

    void getImgData(std::vector<cv::Mat> &B, std::vector<cv::Mat> &L, std::vector<int> &objN);
    void clustDetect(const std::vector<cv::Mat> &B_tmp, const std::vector<cv::Mat> &L_tmp, const std::vector<int> &objN_tmp, const int &fN, const cv::Size &S, const int &t0, const bool &direct, std::vector<cv::Mat> &bw, feature_data &objList, const float &Thr);
    void imgOverlap(const cv::Mat &Bw, cv::Mat &L0, int &cellsIDc, feature_data &featsL1, const float &Thr);
    void clustMasking(const feature_data &objList, const cv::Size &S, const int &fN, std::vector<cv::Mat> &bw);

    void clustDetectRun(const std::vector<cv::Mat> &images, std::vector<cv::Mat> &objMaps, feature_data &objList_final, const float &Thr, const bool &FLAG_DEBUG);
    void assign_ibp_class_to_region(std::vector<std::vector<Region>> &regions, const std::vector<cv::Mat> &objMaps, const int &DELAY, const std::string &INPUTDIR_COLOR, const std::string &OUTDIR, const bool &FLAG_DEBUG);

    /// functions to split the detected clusters based on an watershed segmentation
    int setRandomMarker(const cv::Mat &src, cv::Mat &dst, const int &num_regions);
    int segment_singleCells(const cv::Mat &src_original, const cv::Mat &src_binary, cv::Mat &dst, const int &canny_threshold);
    int distanceTransformation(const cv::Mat &src, cv::Mat &dst);
    void watershed_on_markers(const cv::Mat &src, cv::Mat &markers_pre, cv::Mat &dst);
    void cluster_splitting(std::vector<cv::Mat> &images_gray, std::vector<cv::Mat> &images_binary, feature_data objList, std::vector<cv::Mat> &images_splitted, const int &canny_threshold, const std::string &path, const int &N_THREADS, const bool &FLAG_DEBUG);

} // namespace ibp_cluster_detection
