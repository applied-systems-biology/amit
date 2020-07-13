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

#include "ibp_cluster_detection.h"
#include <iostream>
#include <math.h>  
#include <opencv2/imgproc.hpp>
#include <numeric>
#include <chrono>
#include "IOputs.h"
#include "Outputs.h"
#include "ImageProcessingToolbox.h"
#include "visualize.h"

using std::chrono::system_clock;
using std::chrono::duration;
namespace IPT = ImageProcessingToolbox;


namespace ibp_cluster_detection 
{   
 
    /**
     * Pre-process image data by remove small objects and keep seperated images as labels for all consecutive frames.
     *      
     * @param B all binary images.
     * @param L resulting label images.
     * @param objN resulting number of regions per frame. 
     */
    void getImgData(std::vector<cv::Mat> &B, std::vector<cv::Mat> &L, std::vector<int> &objN){
        
        for (size_t i = 0; i < B.size(); i++) {

            // label all objects and extract the number of objects per frame
            cv::Mat1b tmp_bwlabel;
            int n_labels;
            IPT::bwlabel(B[i], tmp_bwlabel, n_labels, 8);

            L.push_back(tmp_bwlabel);
            
            // subtract 1 because of background will be counted 
            objN.push_back(n_labels - 1);

        }

    }

    /**
     * Perform cluster detection in forward or backwards projection.
     * 
     * @param B_tmp all binary images.
     * @param L_tmp all images within labeled objects.
     * @param objN_tmp all number of regions per frame.
     * @param fN number of frames to analyse.
     * @param S sampling of the input images.
     * @param t0 number of the frame where the detections will be started from (usually =0).
     * @param direct direction of the projection: forward (direct=true) or backwards (direct=false).
     * @param bw resulting image-masks for each object being a cluster. 
     * @param objList object information for all frames.
     * @param Thr threshold for pixel fraction in percents.
     */
    void clustDetect(const std::vector<cv::Mat> &B_tmp, const std::vector<cv::Mat> &L_tmp, const std::vector<int> &objN_tmp, const int &fN, const cv::Size &S, const int &t0, const bool &direct, std::vector<cv::Mat> &bw, feature_data &objList, const float &Thr) {
        
        /// distinguish between forward and backward computation
        std::vector<cv::Mat> B( B_tmp );
        std::vector<cv::Mat> L( L_tmp );
        std::vector<int> objN( objN_tmp );
        
        if (!direct) {
            std::reverse(B.begin(), B.end());
            std::reverse(L.begin(), L.end());
            std::reverse(objN.begin(), objN.end());               
        }           

        ///// cluster detect init

        //  initial labeling, each object on 1st image will get unique ID
        cv::Mat1w L0 = L[t0];
        int num0 = objN[t0];
        
        /// initialisation: list of linear indices of the pixels in the region, returned as a p-element vector
        // about array indexation in mlab run command "doc ind2sub"
        std::vector<std::vector<cv::Point>> objList0_regionProps;
        IPT::regionprops(L0, objList0_regionProps, "PixelIdxList");

        // create and init struct to keep variables for all ROI's
        feature_data objList0( objList0_regionProps );
        // fill with value: {0}
        objList0.ID = std::vector<int> (num0, 0);
        // fill with vector: when a cluster occurs vector has more than one value, default value {0}
        objList0.IDin = std::vector<std::vector<int>> (num0, std::vector<int> {0} );
        objList0.IDout = std::vector<int> (num0, 0);
        // fill with value: {1}
        objList0.In = std::vector<int> (num0, 1);
        objList0.Out = std::vector<int> (num0, 1);
        objList0.num_regions = std::vector<int> (num0, 1);
        // fill with value: {t0}
        objList0.t = std::vector<int> (num0, t0);

        for(int i = 0; i < num0; i++) {

            /// ensure that objects in first frame get unique IDs
            std::vector<int> u, u_pxl_values;
            for (size_t j = 0; j < objList0_regionProps[i].size(); j++) {
                u_pxl_values.push_back( (int) L0.at<ushort>( objList0_regionProps[i][j] ) );
            }

            IPT::tabulate u_tmp( u_pxl_values );

            objList0.ID[i] = u_tmp.value[0];
            objList0.IDout[i] = u_tmp.value[0];

        }

        ///// prealloc memory and copy objList0 to objList, where all variables for all objects will be stored
        objList = objList0;
        
        ///// run cycle (from objects objN[t] to objN[t+1])
        
        // from 2nd frame to frame with index fN
        for (int t = t0+1; t < fN; t++) {

            // analysis of overlapping for two consequent images
            feature_data objListCur;

            ibp_cluster_detection::imgOverlap(B[t], L0, num0, objListCur, Thr);
            
            std::vector<int> tmp = std::vector<int> ( objListCur.ID.size() , t);
            objListCur.t = tmp;

            // concatenate infos from current time-frame with the previous ones
            objList.append( objListCur );

        }

        // a 3D array of binary masks
        ibp_cluster_detection::clustMasking(objList, S, fN, bw);   

        // reverse in case of backward
        if (!direct) {
            std::reverse(bw.begin(), bw.end());            
        }             

    }

    /**
     * Analysis of overlapping for two consequent images.
     * 
     * @param Bw binary image of current frame (CV_8UC1 image, input parameter).
     * @param L0 indexed image from previous frame (mlab: L0 = L2, input parameter).
     * @param cellsIDc counter for IDs for each object.
     * @param featsL1 objects history for current frame (output parameter). 
     * @param Thr threshold for pixel fraction in percents.
     */
    void imgOverlap(const cv::Mat &Bw, cv::Mat &L0, int &cellsIDc, feature_data &featsL1, const float &Thr) {
    
        // labeling of binary image
        cv::Mat L1;
        int num1;
        
        IPT::bwlabel(Bw, L1, num1, 8);
        // subtract 1 because of background will be counted, but discard background
        num1 -= 1;

        // extraction of list of pixels for each object on binary image
        std::vector<std::vector<cv::Point>> featsL1_regionProps;
        IPT::regionprops(L1, featsL1_regionProps, "PixelIdxList");

        featsL1.PixelIdxList = featsL1_regionProps;
        
        // fill with (default) value: {1}
        featsL1.In = std::vector<int> (num1, 1);
        featsL1.Out = std::vector<int> (num1, 1);
        featsL1.t = std::vector<int> (num1, 1);
        featsL1.num_regions = std::vector<int> (num1, 1);
        // fill with vector: when a cluster occurs vector has more than one value
        featsL1.IDin = std::vector<std::vector<int>> (num1, std::vector<int>());
        
        ///// detect in
        /// initialisation of an output indexed image
        cv::Mat L2 = cv::Mat::zeros(Bw.size(), CV_16UC1);

        /// cycle over all objects on binary image Bw
        for (int i = 0; i < num1; i++) {
            /// count unique values of pixels of current object on Bw from L0
            std::vector<int> u, u_pxl_values;
            for (size_t j = 0; j < featsL1.PixelIdxList[i].size(); j++) {
                u_pxl_values.push_back( (int) L0.at<ushort>( featsL1.PixelIdxList[i][j] ) );    
            }

            IPT::tabulate u_tmp( u_pxl_values );

            for(size_t p = 0; p < u_tmp.percent.size(); p++) {
                if(u_tmp.percent[p] > Thr){
                    u.push_back(  u_tmp.value[p] );
                }            
            }

            // if length of vector u
            if ( u.size() == 1 ) 
            {
                // and 1st element of u>0 then cell stay still and ID should not change
                if (u[0] > 0){
                    featsL1.ID.push_back( u[0] );
                    featsL1.IDin[i].push_back( u[0] );
                }
                // and 1st elemen of u==0  then new cell appear
                else {
                    cellsIDc++;
                    featsL1.ID.push_back( cellsIDc );
                    featsL1.IDin[i].push_back( 0 );                                       
                }
            }
            else if ( u.size() == 2 )
            {
                // and 1st element of u==0  then isolated cell move and ID should not change
                if (u[0] == 0){
                    featsL1.ID.push_back( u[1] );
                    featsL1.IDin[i].push_back( u[1] );
                }
                // and 1st element of u>0 then new object move inside cluster
                else{
                    // in that case the ID counter must increment
                    cellsIDc++;
                    // and this object will get new ID
                    featsL1.ID.push_back( cellsIDc );
                    featsL1.In[i] = 2;
                    featsL1.IDin[i].insert( featsL1.IDin[i].end(), u.begin(), u.end() );
                }
            }
            // is > 2 then cells from previous frame formed cluster
            else
            {
                // in that case the ID counter must increment
                cellsIDc++;
                // and this object will get new ID
                featsL1.ID.push_back( cellsIDc );

                // which may contain background pixels
                if (u[0] == 0){
                    // add the number of region when a cluster occur, added by PP
                    featsL1.num_regions[i] = u.size() - 1; // TODO

                    featsL1.In[i] = u.size() - 1;
                    
                    // add all elements from second to last
                    for (size_t k = 1; k < u.size(); k++){
                        featsL1.IDin[i].push_back( u[k] );
                    }
                } 
                // or may not
                else {
                    // add the number of region when a cluster occur, added by PP
                    featsL1.num_regions[i] = u.size(); // TODO

                    featsL1.In[i] = u.size();
                    featsL1.IDin[i].insert( featsL1.IDin[i].end(), u.begin(), u.end() );
                }
        
            } 
            
        }

        
        ///// detect out
        // if cells doesn't interact then output IDs should be the same
        featsL1.IDout = featsL1.ID; 
        
        // however it must be checked with the table of IDs occurence frequency
        IPT::tabulate Tbl_ID( featsL1.ID );

        /// Version 2: looking for IDs which occure more than one time (object from previouse frame splits to few objects in current frame)
        std::vector<int> listID_tmp, listIDin_tmp, featsL1_IDin_tmp, listID;
        
        for(size_t c = 0; c < Tbl_ID.count.size(); c++) {
            if(Tbl_ID.count[c] > 1){
                listID_tmp.push_back( Tbl_ID.value[c] );
            }            
        }
        
        // If there are more that one object have a contribution from one object from previous frame (cluster splitted),
        // then all these objects must receive new IDs

        // iterate through all list of IDin because IDin can contain several values (lists)
        for(size_t obj = 0; obj < featsL1.IDin.size(); obj++) {
            // collect all IDs, where IDin > 1 
            for(size_t id_in = 0; id_in < featsL1.IDin[obj].size(); id_in++) {
                featsL1_IDin_tmp.push_back( featsL1.IDin[obj][id_in]  );
            }
        }

        IPT::tabulate Tbl( featsL1_IDin_tmp );
        
        for(size_t c = 0; c < Tbl.count.size(); c++) {
            if(Tbl.count[c] > 1){
                listIDin_tmp.push_back( Tbl.value[c] );
            }            
        }

        // look for unique IDs by concatenate the two listIDs
        listID = listID_tmp;
        listID.insert( listID.end(), listIDin_tmp.begin(), listIDin_tmp.end() );

        std::sort( listID.begin(), listID.end() );
        listID.erase( std::unique( listID.begin(), listID.end() ), listID.end() );

        // check for Ids > 0, remove if it is not the case
        for(size_t i = 0; i < listID.size(); i++) {
            if(listID[i] <= 0){
                listID.erase( listID.begin() + i );
            }            
        }
        
        // handle for non-empty listID
        if ( ! listID.empty() ){ 
            for(int my_id : listID) {
                // find indices of listID in the ID's - features
                std::vector<int> ids;
                
                for (size_t j = 0; j < featsL1.ID.size(); j++) {
                    if (featsL1.ID[j] == my_id ) {
                        ids.push_back( j );
                    }
                }

                // looking for each of that IDs
                for (size_t j = 0; j < ids.size(); j++) {
                    // and replace them with new IDs
                    cellsIDc++;
                    
                    featsL1.IDout[ ids[j] ] = cellsIDc;
                    featsL1.Out[ ids[j] ] = ids.size();
                }
        
            }
        }

        // an output indexed image for current frame
        for (int i = 0; i < num1; i++){
            for (size_t j = 0; j < featsL1.PixelIdxList[i].size(); j++) {
                L2.at<ushort>( featsL1.PixelIdxList[i][j] ) = featsL1.IDout[i];  
            }            
        }

        // parameter L2 corresponds L0 indexed image from previous frame ( mlab: L0 = L2)
        L0.release();
        L2.copyTo(L0);

    }

    /**
     * Build a 3D array of binary masks where each foreground object indicates a cluster.
     * 
     * @param objList objects history for the current frame.
     * @param S size of image in x and y direction (sampling).
     * @param fN number of frames which corresponds to the 3-dimension like a z-stack.
     * @param bw resulting binary masks which indicate the clusters within (3D mask).  
     */
    void clustMasking(const feature_data &objList, const cv::Size &S, const int &fN, std::vector<cv::Mat> &bw){

        std::vector<int> idxIn;
        
        /// select all objects which have contribution from at least two objects from previous frame
        for (size_t i = 0; i < objList.In.size(); i++){
            if (objList.In[i] > 1) {
                idxIn.push_back ( i );
            }
        }

        // collect all IDs from objects infos
        std::vector<int> idIn;

        for(int elem: idxIn) {
            idIn.push_back( objList.IDout[ elem ] );
        }

        // build a 3D array of binary masks
        bw.clear();
        bw.reserve(fN);
        // init 3D mask with zero images
        for (int i = 0; i < fN; i++){
            bw.push_back( cv::Mat1b::zeros( S ) );
        }

        // concatenate indices
        std::vector<int> ids = idIn;
        
        // fill the 3D-mask
        for (size_t i = 0; i < ids.size(); i++)
        {
            std::vector<int> idx;
            for (size_t j = 0; j < objList.ID.size(); j++)
            {        
                if (objList.IDout[j] == ids[i]){
                    idx.push_back( j );
                }
            }

            for (size_t j = 0; j < idx.size(); j++)
            {
                // select objList[idx[j]].t channel / stack where a cluster occurs at timepoint t
                cv::Mat tmp = bw[objList.t[ idx[j] ]];

                // set corresponding pixel-values where a cluster is located to 1 (255)
                for(cv::Point pixel: objList.PixelIdxList[ idx[j] ]) {
                    tmp.at<uchar>( pixel ) = 255;
                }

                // assign the corresponding stack / mask
                bw[objList.t[ idx[j] ]] = tmp;
            }
            
        }

    }

    /**
     * Perform the main cluster detection method for all given frames with binary images. 
     *
     * @param images binary input images. 
     * @param objMaps resulting masks for all frames where objects with intensity 127 indicates single cells and intensity value 255 to be a cluster (0 = noise).
     * @param objList_final combined list with history information for each object.
     * @param Thr threshold for pixel fraction in percents.
     * @param FLAG_DEBUG verbose parameter.
     */    
    void clustDetectRun(const std::vector<cv::Mat> &images, std::vector<cv::Mat> &objMaps, feature_data &objList_final, const float &Thr, const bool &FLAG_DEBUG){

        // number of frames to analyse
        const int fN = images.size();
        
        // the resolution of input images
        cv::Size S = images[0].size();  
        
        /// pre-processing: load binary images,remove small objects and label the rest of it
        std::vector<cv::Mat> B( images );
        std::vector<cv::Mat> L;

        // number ob objects per each frame
        std::vector<int> objN;         
        
        ibp_cluster_detection::getImgData(B, L, objN);
        
        // start time measurement
        auto start = system_clock::now();
        
        ///// cluster detection: forward passage, from the 1st frame to the last one
        feature_data objList;

        std::cout << " [3.0.1] clustDetectRun: start forward cluster detection" << std::endl;

        std::vector<cv::Mat> bwF;
        ibp_cluster_detection::clustDetect(B, L, objN, fN, S, 0, true, bwF, objList, Thr);

        std::cout << " [3.0.2] clustDetectRun: start backward cluster detection" << std::endl;
        
        /// cluster detection: backward passage, from the last frame to the 1st one
        feature_data objList_backwards;

        std::vector<cv::Mat> bwB;
        ibp_cluster_detection::clustDetect(B, L, objN, fN, S, 0, false, bwB, objList_backwards, Thr);

        std::cout << " [3.0.3] clustDetectRun: finish cluster detection - perform overlay" <<  std::endl;
        std::cout << " [3.0.4] clustDetectRun: bwF.size: " << bwF.size() << "\tbwB.size: " << bwB.size() << std::endl;

        ///// clear cv::Mat vectors that are no longer needed to free memory  
        L.clear();

        // array of masks of clusters for each frame = element by element logical multiplication 
        // of two masks for each time frame
        std::vector<cv::Mat> clustM;
        clustM.reserve( fN );

        for (size_t i = 0; i < bwF.size(); i++){
            cv::Mat dst;

            cv::bitwise_or(bwF[i], bwB[i], dst);
            clustM.push_back( dst );
        }

        ///// clear cv::Mat vectors that are no longer needed to free memory 
        bwF.clear();
        bwB.clear();
        
        // array of masks of single objects for each frame = element by element logical multiplication of
        // all object masks and inverted masks of clusters for each time frame
        std::vector<cv::Mat> singlM;
        singlM.reserve( fN );

        for (size_t i = 0; i < B.size(); i++){
            cv::Mat clustM_i_inverted, dst;
            cv::bitwise_not( clustM[i], clustM_i_inverted );

            cv::bitwise_and(B[i], clustM_i_inverted, dst);
            singlM.push_back( dst );            
        }

        ///// clear cv::Mat vectors that are no longer needed to free memory  
        B.clear();        

        // array of objects maps for each frame: 0 - background, 127 - single, 255 - cluster
        // uint8 - conversion to 8 bit integer type: A.*B - element-by-element multiplication
        std::vector<cv::Mat> objMap;
        objMap.reserve( fN );

        for (size_t i = 0; i < singlM.size(); i++){
            cv::Mat tmp, dst;
            
            tmp = singlM[i] - cv::Scalar(128) ; 

            cv::bitwise_or(clustM[i], tmp, dst);

            objMap.push_back( dst );
        }

        ///// fill info about each object
        std::cout << " [3.0.5] clustDetectRun: objMap.size(): " <<  objMap.size() << "\tID.size(): " <<  objList.ID.size() << std::endl;
                
        for (size_t i = 0; i < objList.ID.size(); i++){

            cv::Mat objMapT = objMap[ objList.t[i] ];
             
            // assign type of region by assign one (first) pixel-value of corresponding pixel             
            int pixel_value = (int) objMapT.at<uchar>( objList.PixelIdxList[i][0] );

            // change pixel value in range from {127,255} to {1,2} = { singlecells, cluster}
            if(pixel_value == 127)
                pixel_value = 1;
            else 
                pixel_value = 2;
            
            objList.type.push_back( pixel_value );
       
        }

        ///// determine the number of regions per cluster (single cell = 1)
        objList.getNumberRegions( objList_backwards );
       
        // print total calculation time
        auto end = system_clock::now();
        const double elapsed_seconds = duration<double>(end-start).count();
        std::cout << " [3.0.6] clustDetectRun: end, took\t" << elapsed_seconds << " seconds" << std::endl;

        // copy variables for output
        objMaps = objMap;
        objList_final = objList;

        if (FLAG_DEBUG){
            // print out the object list for debugging purpose ...
            std::cout << "objList.ID.size(): " << objList.ID.size() << std::endl;
            std::cout << "\ni\t;PixelIdxList\t; ID\t; IDin\t; IDout\t; In\t; Out\t; t\t; type\t; num_regions" << std::endl;
            std::cout << "__________________________________________________________________" << std::endl;
            for(size_t idx = 0; idx < objList.ID.size(); idx++) {
                std::cout << idx << "\t; " << objList.PixelIdxList[idx].size() << "\t\t; " << objList.ID[idx] << "\t; ";
                for(int id_in_tmp : objList.IDin[idx]) {
                    std::cout << "[" << id_in_tmp << "]";
                }
                std::cout << "\t; " << objList.IDout[idx] << "\t; " << objList.In[idx] << "\t; "
                          << objList.Out[idx] << "\t; " << objList.t[idx]+1 << "\t; " << objList.type[idx] << "\t; "
                          << objList.num_regions[idx] << std::endl;
            }
            
        }        

    }

    /**
     * Assign the class for each ROI, where the information is taken from the cluster-detection-masks.
     * 
     * @param regions regions-per-frames, which contain the variable "klass", that will be assigned to noise, single cells or cluster.
     * @param objList image masks for all frames where objects with intensity 127 indicates single cells and intensity value 255 to be a cluster (0 = noise).
     * @param DELAY optional delay between binary and original images.
     * @param INPUTDIR_COLOR path where the original gray-scaled images are stored.
     * @param OUTDIR path where the final images with their corresponding classified objects will be stored.
     * @param FLAG_DEBUG verbose parameter.
     */ 
    void assign_ibp_class_to_region(std::vector<std::vector<Region>> &regions, const std::vector<cv::Mat> &objMaps, const int &DELAY, const std::string &INPUTDIR_COLOR, const std::string &OUTDIR, const bool &FLAG_DEBUG){

        // iterate about all frames
        for (size_t i = 0; i < regions.size(); i++){
            
            if (FLAG_DEBUG){
                std::cout << " #ROI's - frame: "  << regions[i].size() << " / " << i << std::endl;    
            }            
            
            std::vector<Region>::iterator reg;
            // iterate about all regions            
            for(reg = regions[i].begin(); reg != regions[i].end(); reg++){

                // select the first pixel of all pixel within the region (center is sometimes not within the region itself)
                cv::Vec2i region_point = reg->region_pixels[0];
                // objMaps[i] works, because of the same counter-sequence of objects: column-wise
                int region_value = objMaps[i].at<uchar>( region_point );

                // assign center_value of region accordingly to the class of regions: 0 (should not appear) = background; 127 = single cell; 255 = cluster
                if (region_value == 0){
                    // noise / background   = 0
                    reg->setClass( cell_class::immune::NOISE );
                }
                else if(region_value == 127){
                    // single cells         = 1  
                    reg->setClass( cell_class::immune::SINGLE );
                }
                else if(region_value == 255){
                    // cluster              = 2
                    reg->setClass( cell_class::immune::CLUSTER );
                }
                else{
                    // ERROR: default assignment => default assignment of class <single cell>
                    std::cout << " CAUTION: center-value of detected ROI is neither associated to background,singlecell,cluster = (0,127,255)" << std::endl;
                    reg->setClass( cell_class::immune::SINGLE );
                }
                
                if (FLAG_DEBUG){
                    cell_class::immune myclass = reg->getClass();
                    std::cout << "\tregion-point: " << region_point << "\tregion-value: " << region_value << "\tclass of ROI: " << cell_class::getCellClass(myclass) << std::endl;
                }

            }
            
        }

        // save images with drawn outline dependent on their class 
        if(FLAG_DEBUG){
			std::string file1 = OUTDIR + "/3_ibp_region_classification/";
            
            outputs::showClassifiedSegmentedRegions2D(file1, INPUTDIR_COLOR , DELAY, regions, true, false);
		}

    }

    /**
     * Determine a specified number of (random) seed points within the cluster.
     * 
     * For 1-5 values the following seed points will be generated: top-most, bottom-most, left-most, right-most, center point.
     * If there are more than 5 points the rest of the points will be generated randomly within the cluster. 
     *
     * @param src binary input cut-out image of the cluster.
     * @param dst resulting binary cut-out image within the generated seed points.
     * @param num_regions of seed points which should be generated. 
     * @return number of distinct regions after generating (random) seed points.
     */
    int setRandomMarker(const cv::Mat &src, cv::Mat &dst, const int &num_regions){

        dst = cv::Mat::zeros(src.size(), CV_8UC1);

        /// erode the original cluster ROI to ensure creating random seed points inside the original ROI
        cv::Mat src_erode;
        cv::Mat kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(3,3));
        cv::erode(src, src_erode, kernel);

        std::vector<std::vector<cv::Point> > contours;
        std::vector<cv::Vec4i> hierarchy;
        cv::findContours(src_erode, contours, hierarchy, cv::RetrievalModes::RETR_LIST, cv::CHAIN_APPROX_SIMPLE);

        if ( contours.size() == 0 ) {
            std::cout << " WARNING: exception-handling: no contoures inside function <setRandomMarker> , return 1" << std::endl;
            return 1;
        }

        std::vector<cv::Point> cnt = contours[0];

        /// calculate the extreme 4-extreme points: [ left, right, top, bottom ]
        auto left_right = std::minmax_element(cnt.begin(), cnt.end(), [](cv::Point const& a, cv::Point const& b){
            return a.x < b.x;
        });

        auto top_bottom = std::minmax_element(cnt.begin(), cnt.end(), [](cv::Point const& a, cv::Point const& b){
            return a.y < b.y;
        });

        /// CAUTION: always go one pixel into the cluster to ensure be inside the ROI
        cv::Point leftMost = cv::Point(left_right.first->x+1 , left_right.first->y );
        cv::Point rightMost = cv::Point(left_right.second->x-1 , left_right.second->y );
        cv::Point topMost = cv::Point(top_bottom.first->x , top_bottom.first->y+1 );
        cv::Point bottomMost = cv::Point(top_bottom.second->x , top_bottom.second->y-1 ) ;

        /// calculate if the cluster is wider than higher
        int width = abs( left_right.first->x - left_right.second->x );
        int height = abs( top_bottom.first->y - top_bottom.second->y );

        cv::Rect rect(cv::Point(), src.size());

        std::vector<cv::Point> seed_points;

        /// EXCEPTION case when there is just on distinct is detected (if clear border switched of)
        if ( num_regions == 1 ){
            /// assign centroid of cluster
            cv::Moments m = moments(src,true);
            cv::Point center(m.m10/m.m00, m.m01/m.m00);

            /// take care that this center-point it inside the cluster-region
            if (rect.contains(center) && src.at<uchar>( center ) == 0) {
                /// select the top point - 1 in y direction when the center is outside the cluster ROI
                seed_points.push_back( topMost );
            } else {
                /// center is inside the cluster ROI
                seed_points.push_back( center );
            }

        }
        /// distinguish for random seed points ( cluster is wider than high )
        else if ( width > height){
            seed_points.push_back( leftMost );
            seed_points.push_back( rightMost );

            /// collect for 3 regions additionally top most point
            if (num_regions == 3){
                seed_points.push_back( topMost );
            }
            /// collect for more than 3 regions additionally most top and bottom point
            else if(num_regions > 3) {
                seed_points.push_back( topMost );
                seed_points.push_back( bottomMost );
            }

        }
        ///  cluster is higher than wide
        else {
            /// collect most top and bottom point
            seed_points.push_back( topMost );
            seed_points.push_back( bottomMost );

            /// collect for 3 regions additionally left most point
            if (num_regions == 3){
                seed_points.push_back( leftMost );
            }
            /// collect for more than 3 regions additionally most left and right point
            else if(num_regions > 3) {
                seed_points.push_back( leftMost );
                seed_points.push_back( rightMost );
            }
        }

        /// add center of cluster for more than 4 needed seed points
        if (num_regions > 4){

            cv::Mat src_center;
            src.copyTo(src_center);

            cv::Point center;

            /// take care that this center-point it inside the cluster-region, erode as long as center is inside the ROI
            while (rect.contains(center) && src_erode.at<uchar>( center ) == 0) {
                cv::erode(src_center, src_center, kernel);

                cv::Moments m = moments(src_center,true);
                center = cv::Point(m.m10/m.m00, m.m01/m.m00);
            }

            seed_points.push_back( center );

            /// more than 5 objects: collect all pixels within the object
            if (num_regions > 5){

                /// select random seed on eroded cluster ROI to ensure creating random seed points inside the original ROI
                cv::Mat labels, stats, centroids;
                (void)connectedComponentsWithStats(src_erode, labels, stats, centroids);
                
                std::vector<cv::Point> all_pixels;

                /// just select label = 1 (must be just one label because of just cluster is included and 0 = background)
                int i = 1;

                /// get bounding rect
                int left =  stats.at<int>(i, cv::CC_STAT_LEFT) ;
                int top = stats.at<int>(i, cv::CC_STAT_TOP);
                int width = stats.at<int>(i, cv::CC_STAT_WIDTH);
                int height = stats.at<int>(i, cv::CC_STAT_HEIGHT);

                int x_end = left + width;
                int y_end = top + height;
                for (int x = left; x < x_end; x++) {
                    for (int y = top; y < y_end; y++) {
                        cv::Point p(x, y);

                        if ( i == labels.at<int>(p) ) {
                            all_pixels.push_back(p);
                        }
                    }

                }

                /// after that add as much as possible randomly
                while ( (int)seed_points.size() < num_regions ) {
                    /// create random index for extract one pixel of whole cluster
                    int randomIndex = rand() % all_pixels.size();

                    /// check if point already exist is list or is on the contour
                    std::vector<cv::Point>::iterator it = std::find(seed_points.begin(), seed_points.end(), all_pixels[randomIndex] );

                    /// take also care that random seed point is inside the cluster and not already exist
                    if ( it != seed_points.end() or (rect.contains(all_pixels[randomIndex]) && src_erode.at<uchar>( all_pixels[randomIndex] ) == 0) )
                        continue;
                    else
                        seed_points.push_back( all_pixels[randomIndex] );

                }

            }

        }

        for (const auto & seed_point : seed_points) {
            cv::circle(dst, seed_point, 1, cv::Scalar(255), 1);
        }

        cv::Mat labels;
        int num_distinct_regions = cv::connectedComponents(dst, labels, 4) - 1;

        return num_distinct_regions;

    }

    /**
     * Segment the inner cells (within a cluster) of the cut-out image with the original ibp-segmentation.
     *
     * @param src_original gray-scaled input image.
     * @param src_binary already segmented cluster of the cut-out image.
     * @param dst resulting binary cut-out image.
     * @param canny_threshold upper threshold for canny edge detector, the higher the more sensitive the segmentation is.
     * @return number of distinct regions after segmentation.
     */
    int segment_singleCells(const cv::Mat &src_original, const cv::Mat &src_binary, cv::Mat &dst, const int &canny_threshold){

        cv::Mat adapthisteq, blur, canny, bw_threshold, bw_closeFill, bw_overlay, bwareaopen;

        /// contrast enhancement via uniform distribution
        IPT::adapthisteq(src_original, adapthisteq);

        /// blur image to reduce noise
        cv::blur( adapthisteq, blur, cv::Size(3,3) );

        /// detect edges using canny: threshold was figured out via visual::fincContours
        std::vector<std::vector<cv::Point> > contours_canny;
        std::vector<cv::Vec4i> hierarchy_canny;

        cv::Canny(blur , canny, canny_threshold, 255, 3 );
        /// find contours
        cv::findContours(canny, contours_canny, hierarchy_canny, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );

        /// draw contours
        bw_threshold = cv::Mat::zeros(canny.size(), CV_8UC1 );
        for(size_t i = 0; i< contours_canny.size(); ++i ) {
            cv::drawContours(bw_threshold, contours_canny, i, cv::Scalar(255), 2, 8, hierarchy_canny, 0, cv::Point() );
        }

        /// normalize all values in range [ 0 , 255 ]
        cv::normalize(bw_threshold, bw_threshold, 0, 255.0, cv::NORM_MINMAX);

        /// perform morphology - closing for each connected component (cell) separate
        cv::Mat1b bw_label_tmp;
        int n_labels;

        IPT::bwlabel(bw_threshold, bw_label_tmp, n_labels, 4);

        std::vector<float> label_values = IPT::unique(bw_label_tmp, true);

        cv::Mat bw_label;
        bw_label_tmp.copyTo(bw_label);
        bw_label.convertTo(bw_label, CV_8UC1);

        cv::Mat combined_components = cv::Mat::zeros(bw_label.size(), CV_8UC1);

        /// iterate over all single component / label
        for (float label_value : label_values) {

            /// skip background
            if ((int)label_value == 0) {
                continue;
            }

            cv::Mat single_component = cv::Mat::zeros(bw_label.size(), CV_8UC1);

            /// iterate over the label image
            for (int i = 0; i < bw_label.rows; i++) {
                for (int j = 0; j < bw_label.cols; j++) {

                    int editValue = bw_label.at<uchar>(i, j);

                    if ( editValue == (int)label_value ) {
                        single_component.at<uchar>(i,j) = 255;
                    }
                }
            }

            /// perform the closing on the single ROI
            cv::Mat kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(5,5));

            cv::Mat single_closing;
            cv::morphologyEx(single_component, single_closing, cv::MORPH_CLOSE, kernel, cv::Point(-1, -1), 1);

            /// fill the remaining holes of single component
            cv::Mat single_fill;
            IPT::imfill(single_closing, single_fill);

            /// remove objects smaller than x pixels
            cv::Mat single_areaopen;
            IPT::bwareaopen(single_fill, single_areaopen, 250);

            /// combine single component together
            cv::bitwise_or(combined_components, single_areaopen, combined_components);

        }

        /// binary bw_overlay with original cluster segmentation to avoid "over-segmentation" and keep into the segmented region
        cv::bitwise_and(src_binary, combined_components, dst);

        /// compute distinct regions (-1 because the background will be counted)
        cv::Mat labels;
        int num_distinct_regions = cv::connectedComponents(dst, labels, 4) - 1;

        return num_distinct_regions;

    }

    /**
     * Calculate the distance transormation for each foreground pixel and it's distance to the next background pixel.
     *
     * @param src binary input cut-out image.
     * @param dst resulting binary cut-out image with distance thresholding.
     * @return number of distinct regions after distance transformation.
     */
    int distanceTransformation(const cv::Mat &src, cv::Mat &dst){

        /// perform the distance transform algorithm
        cv::Mat dist, threshold, opening, labels;

        cv::distanceTransform(src, dist, cv::DIST_WELSCH, cv::DIST_MASK_PRECISE);

        /// normalize the distance image for range = {0.0, 1.0} so we can visualize and threshold it
        cv::normalize(dist, dist, 0, 1.0, cv::NORM_MINMAX);

        /// threshold to obtain the peaks, this will be the markers for the foreground objects, convert back to range [ 0,255 ]
        cv::threshold(dist, threshold, 0.4, 1.0, cv::THRESH_BINARY);
        threshold.convertTo(threshold, CV_8UC1, 255, 0);

        /// opening to remove edge-artifacts as well as to small regions
        cv::Mat kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(3,3));

        cv::morphologyEx(threshold, opening, cv::MORPH_OPEN, kernel);

        opening.copyTo(dst);

        int num_labels = connectedComponents(dst, labels, 4) - 1;

        return num_labels;

    }

    /**
     * Perform the watershed segmentation on the provided seed regions inside the cluster.
     *
     * @param src binary input cut-out image.
     * @param markers_pre binary cut-out image with corresponding seed points/markers for the watershed segmenation.
     * @param dst resulting cut-out image after watershed segmentation.
     */
    void watershed_on_markers(const cv::Mat &src, cv::Mat &markers_pre, cv::Mat &dst){

        /// convert binary and final image to "BGR"
        cv::Mat imgResult;
        src.copyTo(imgResult);
        cv::cvtColor(imgResult, imgResult, cv::COLOR_GRAY2BGR);

        /// create the CV_8U version of the distance image, it is needed for findContours()
        cv::Mat dist_8u;
        markers_pre.convertTo(dist_8u, CV_8U);
        /// find total markers in distanceTransformed image
        std::vector<std::vector<cv::Point>> contours;

        cv::findContours(dist_8u, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

        /// create the marker image for the watershed algorithm
        cv::Mat markers = cv::Mat::zeros(markers_pre.size(), CV_32S);
        /// draw the foreground markers
        for (size_t i = 0; i < contours.size(); i++){
            cv::drawContours(markers, contours, static_cast<int>(i), cv::Scalar(static_cast<int>(i)+1), -1);
        }

        /// draw the background marker
        cv::circle(markers, cv::Point(5,5), 3, cv::Scalar(255), -1);

        /// perform the watershed algorithm
        cv::watershed(imgResult, markers);

        cv::Mat mark;
        markers.convertTo(mark, CV_8U);
        cv::bitwise_not(mark, mark);

        /// image looks like at that point, generate random colors
        std::vector<cv::Vec3b> colors;
        for (size_t i = 0; i < contours.size(); i++){
            int b = cv::theRNG().uniform(0, 256);
            int g = cv::theRNG().uniform(0, 256);
            int r = cv::theRNG().uniform(0, 256);
            colors.push_back(cv::Vec3b((uchar)b, (uchar)g, (uchar)r));
        }

        /// create the result image and fill labeled objects with random colors
        dst = cv::Mat::zeros(markers.size(), CV_8UC3);

        for (int i = 0; i < markers.rows; i++)
        {
            for (int j = 0; j < markers.cols; j++)
            {
                int index = markers.at<int>(i,j);
                if (index > 0 && index <= static_cast<int>(contours.size()))
                {
                    dst.at<cv::Vec3b>(i,j) = colors[index-1];
                }
            }
        }

        /// convert back to gray-scaled image and threshold it
        cv::cvtColor(dst, dst, cv::COLOR_BGR2GRAY);

        cv::threshold(dst, dst, 1.0, 255, cv::THRESH_BINARY);

    }

    /**
     * Split the detected clusters based on an watershed segmentation.
     * 
     * Identify the detected clusters position and perform the watershed on a cut-out subwindow on such cluster.
     * (0) Perform contrast enhancement and noise removal as used for the ibp::segmentation
     * (1) Extract cluster from ibp_cluster_detection
     * (2) Create a mask image and build overlay, that the later bounding box will be taken just on the cluster ROIs
     * (3) Extract and align bounding-box of gray and binary image for each cluster
     * (4) Ibp segment inner cell to get number of markers
     * (5) Extract marker by distance transformation as long as there are so many distinct regions as known for the corresponding cluster
     * (6) Check whether there a still not the correct number of marker, if yes select randomly based on cluster size
     * (7) Split the cluster with the watershed segmentation algorithm
     * (8) Subtract image (1) from the original image and insert the splitted cluster image from (2) to get the final image
     *
     * @param images_gray original gray images
     * @param images_binary original binary images
     * @param objList list with all information about each region based on the ibp cluster detection
     * @param images_dst binary images with splitted cluster
     * @param canny_threshold threshold to detect inner cell membrane within a detected cluster, the higher the less contours remains
     */
    void cluster_splitting(std::vector<cv::Mat> &images_gray, std::vector<cv::Mat> &images_binary, feature_data objList, std::vector<cv::Mat> &images_dst, const int &canny_threshold, const std::string &path, const int &N_THREADS, const bool &FLAG_DEBUG){

        /// create directory to store images
        if(FLAG_DEBUG){
            (void) io::create_directory(path);
        }

        ///// iterate over all images - multi-threaded processing/////
        #pragma omp parallel for firstprivate(images_binary) num_threads(N_THREADS)
        for (size_t t = 0; t < images_binary.size(); ++t) {

            /// (0) Perform contrast enhancement and noise removal as used for the ibp::segmentation
            cv::Mat sd1 = IPT::createKernel("disk3");
            cv::Mat sd2 = IPT::createKernel("disk5");

            cv::Mat Io, Jt, J1, Ic, Jb, img_gray;
            /// sharpen and background correction step 1 (opening,subtract,add)
            IPT::imopen(images_gray[t], Io, sd2);
            IPT::imsubtract(images_gray[t], Io, Jt);
            IPT::imadd(images_gray[t], Jt, J1);

            /// sharpening and background correction step 2 (close, 2xsubtract)
            IPT::imclose(images_gray[t], Ic, sd1);
            IPT::imsubtract(Ic, images_gray[t], Jb);
            IPT::imsubtract(J1, Jb, img_gray);

            cv::Mat img_dst;
            images_binary[t].copyTo(img_dst);

            /// (1) Extract cluster from ibp_cluster_detection
            std::vector<int> idx_cluster;
            std::vector<int> num_regions;

            objList.getClusterIndex(t, idx_cluster, num_regions);

            /// DEBUG : plot the resulting cluster split on the original grayscaled images
            /// create an colored debug image with drawn contours for clusters (red) and single cells (blue)
            cv::Mat img_debug;
            if (FLAG_DEBUG){
                images_gray[t].copyTo(img_debug);

                std::vector<std::vector<cv::Point>> contours_all_cells;
                cv::findContours(img_dst, contours_all_cells, cv::RETR_LIST, cv::CHAIN_APPROX_SIMPLE);

                cv::cvtColor(img_debug, img_debug, cv::COLOR_GRAY2BGR);
                /// draw inner of the original cluster
                for (size_t c_img_debug = 0; c_img_debug < contours_all_cells.size(); c_img_debug++)
                    cv::drawContours(img_debug, contours_all_cells, c_img_debug, cv::Scalar(255, 0, 0), 1);

            }            

            ///// iterate over all cluster and draw the regions /////
            for (size_t idx = 0; idx < idx_cluster.size(); ++idx) {

                int c = idx_cluster[idx];

                /// (2) Create a cluster-image-mask and build overlay, that the later bounding box will be taken just on the one cluster ROI
                cv::Mat img_mask = cv::Mat::zeros(images_binary[t].size(), CV_8UC1);
                cv::drawContours(img_mask, objList.PixelIdxList, c, cv::Scalar(255), cv::FILLED);

                cv::Mat img_binary( images_binary[t].size(), CV_8UC1, cv::Scalar(0) );
                cv::bitwise_and(images_binary[t], img_mask, img_binary);

                /// (3) Extract and align bounding-box of gray and binary image for each cluster
                cv::Rect myrect = boundingRect(objList.PixelIdxList[c]);

                cv::Mat ROI_gray = img_gray(myrect);
                cv::Mat ROI_binary = img_binary(myrect);

                /// padding surround the window to prevent artifacts in watershed-segmentation and inner cell segmentation, later cutout again
                int padding = 4;

                /// get the value in the top-left corner for padding constant value
                int pad_value = (int) ROI_gray.at<ushort>(cv::Point(0,0));
                cv::Mat ROI_gray_pad;
                cv::copyMakeBorder(ROI_gray, ROI_gray_pad, padding, padding, padding, padding, cv::BORDER_CONSTANT, cv::Scalar(pad_value));

                cv::Mat ROI_binary_pad;
                cv::copyMakeBorder(ROI_binary, ROI_binary_pad, padding, padding, padding, padding, cv::BORDER_CONSTANT, cv::Scalar(0));

                /// prepare second ROI_marker directly for the loop
                cv::Mat ROI_markers, ROI_markers_dT;
                ROI_binary_pad.copyTo(ROI_markers_dT);

                int num_distinct_regions;
                int num_distinct_regions_dT = -1;

                /// (4) Ibp segment inner cell to get number of markers                 
                num_distinct_regions = ibp_cluster_detection::segment_singleCells(ROI_gray_pad, ROI_binary_pad, ROI_markers, canny_threshold);

                /// CAUTION: exception handling when number of regions inside cluster = 0 => assign number of already detected markers
                if ( num_regions[idx] == 0){
                    num_regions[idx] = num_distinct_regions;
                }

                int ITERATIONS = 0;

                if (FLAG_DEBUG){
                    std::cout << "\tINFO: # ROI's - marker [ actual / target ]: [ " << num_distinct_regions << " / " << num_regions[idx] << " ] after segment_singleCells [frame]: " << t << std::endl;
                }

                /// (5) Extract marker by distance transformation as long as there are so many distinct regions as known for the corresponding cluster
                while ( num_distinct_regions != num_regions[idx] or num_distinct_regions_dT != num_regions[idx] ) {

                    /// continue with markers from single-cell-segmentation, abort if needed number is reached
                    num_distinct_regions = ibp_cluster_detection::distanceTransformation(ROI_markers, ROI_markers);
                    if (FLAG_DEBUG) {
                        std::cout << "\tINFO: # ROI's - marker    [ actual / target ]: [ " << num_distinct_regions << " / " << num_regions[idx] << " ]\t => distanceTransformation [frame]: " << t << std::endl;
                    }
                    if (num_distinct_regions == num_regions[idx]) break;

                    /// create markers from original binary cluster ROI, abort if needed number is reached
                    num_distinct_regions_dT = ibp_cluster_detection::distanceTransformation(ROI_markers_dT, ROI_markers_dT);
                    if (FLAG_DEBUG) {
                        std::cout << "\tINFO: # ROI's - marker_dt [ actual / target ]: [ " << num_distinct_regions_dT << " / " << num_regions[idx] << " ]\t => distanceTransformation [frame]: " << t << std::endl;
                    }
                    if (num_distinct_regions_dT == num_regions[idx]) break;

                    /// ABORT criterion
                    if (ITERATIONS > 10) {
                        if (FLAG_DEBUG){
                            std::cout << "\tABORT: marker-thinning\t[frame]: " << t << std::endl;
                        }
                        break;
                    }

                    ITERATIONS++;
                }

                /// select the marker ROI with the correct number of distinct region if available
                if ( num_distinct_regions != num_regions[idx] and num_distinct_regions_dT == num_regions[idx] ) {
                    ROI_markers = ROI_markers_dT;
                }

                /// (6) Check whether there a still not the correct number of marker, if yes select randomly based on cluster size
                if ( num_distinct_regions != num_regions[idx] and num_distinct_regions_dT != num_regions[idx] ){
                    
                    if (FLAG_DEBUG){
                        std::cout << "\tINFO: # ROI's - random-marker [ actual / target ]: [ " << num_distinct_regions << " / " << num_regions[idx] << " ]\t => before setRandomMarker [frame]: " << t << std::endl;
                    }

                    try {
                        num_distinct_regions = ibp_cluster_detection::setRandomMarker(ROI_binary_pad, ROI_markers, num_regions[idx]);
                    } catch (const std::exception& e) {
                        std::cout << e.what() << std::endl;
                        num_distinct_regions = 1;
                    }

                    if (FLAG_DEBUG){
                        std::cout << "\tINFO: # ROI's - random-marker [ actual / target ]: [ " << num_distinct_regions << " / " << num_regions[idx] << " ]\t => after setRandomMarker [frame]: " << t << std::endl;
                    }
                    
                }

                /// perform the watershed segmentation on the final markers
                cv::Mat ROI_splitted, ROI_dst;

                /// (7) Split the cluster with the watershed segmentation algorithm
                ibp_cluster_detection::watershed_on_markers(ROI_binary_pad, ROI_markers, ROI_splitted);

                /// cut out the ROI_binary_pad without the padding and set it back in the original binary image
                cv::Rect myrect_dst = cv::Rect(0 + padding, 0 + padding, myrect.width , myrect.height);
                ROI_dst = ROI_splitted(myrect_dst);

                /// extract just unwanted neighbour cell(s)
                cv::Mat ROI_binary_justNeighbours;
                cv::bitwise_xor(ROI_binary, ROI_binary, ROI_binary_justNeighbours);

                /// overlay of unwanted neighbour cell(s) and split cluster
                cv::Mat ROI_dst_merged;
                cv::bitwise_or(ROI_dst, ROI_binary_justNeighbours, ROI_dst_merged);

                /// (8) Subtract image (1) from the original image and insert the splitted cluster image to get the final image
                ROI_dst_merged.copyTo( img_dst(myrect) );

                if (FLAG_DEBUG){
                    cv::Mat ROI_drawn_outline, ROI_debug_merged;
                    ROI_gray_pad.copyTo(ROI_drawn_outline);
                    cv::cvtColor(ROI_drawn_outline, ROI_drawn_outline, cv::COLOR_GRAY2BGR);

                    /// get contour surround the split regions of cluster
                    std::vector<std::vector<cv::Point>> contours_ROI_splitted;
                    cv::findContours(ROI_splitted, contours_ROI_splitted, cv::RETR_LIST, cv::CHAIN_APPROX_SIMPLE);

                    /// draw contour of the split regions in the original gray-scaled ROI (red)
                    for (size_t c_debug_split = 0; c_debug_split < contours_ROI_splitted.size(); c_debug_split++)
                        cv::drawContours(ROI_drawn_outline, contours_ROI_splitted, c_debug_split, cv::Scalar(0, 0, 255), 1);

                    /// get contour surround the neighbours regions (if available) and draw it into the padding image
                    cv::Mat ROI_binary_justNeighbours_pad;
                    ROI_binary_justNeighbours.copyTo(ROI_binary_justNeighbours_pad);
                    cv::copyMakeBorder(ROI_binary_justNeighbours, ROI_binary_justNeighbours_pad, padding, padding, padding, padding, cv::BORDER_CONSTANT, cv::Scalar_(0));

                    std::vector<std::vector<cv::Point>> contours_neighbours;
                    cv::findContours(ROI_binary_justNeighbours_pad, contours_neighbours, cv::RETR_LIST, cv::CHAIN_APPROX_SIMPLE);

                    /// draw contour of the neighbours regions in the original gray-scaled ROI (blue)
                    // for (int c_debug_split = 0; c_debug_split < contours_neighbours.size(); c_debug_split++)
                    // cv::drawContours(ROI_drawn_outline, contours_neighbours, c_debug_split, cv::Scalar(255, 0, 0), 1);

                    /// cut out the split cluster and place it (with original other cells into the box) back into the binary image
                    cv::Mat ROI_debug = ROI_drawn_outline(myrect_dst);

                    /// CAUTION draw lines-contour inside original gray-scaled image
                    ROI_debug.copyTo( img_debug(myrect) );

                }

            }

            images_dst.push_back(img_dst);

            /// store images if debug mode is activate
            if(FLAG_DEBUG){
                io::write_image(path, img_debug, t);
            }

        }

    }


} // namespace ibp_cluster_detection