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

#include "tracking_region_association.h"
#include "Classification.h"
#include "Outputs.h"


namespace tracking_region_association
{

    void setFlatStatus(std::vector<Region> &regions, const bool stat){

        for(size_t i = 0; i < regions.size(); i++){
            regions.at(i).setFlat(stat);
        }
    
    }

    void setEllipsod(std::vector<Region> &regions){
        cv::Size2f s;
        
        for(size_t i = 0; i < regions.size(); i++){
            regions.at(i).getPrincipalAxesMinEnclosingEllipse(s);

            if(s.height < s.width )
                regions.at(i).setMinMajAxe(s.height, s.width);
            else
                regions.at(i).setMinMajAxe(s.width, s.height);
            
        }

    }

    /**
     * set klass - variable of each region
     */
    void assign_gmm_class_to_region(std::vector<std::vector<Region>> &regions, const cv::Mat &labels,  const bool &FLAG_DEBUG){

        int labels_i_w = 0;
        
        for(std::vector<std::vector<Region>>::iterator it = regions.begin(); it != regions.end(); it++){
            
            for(std::vector<Region>::iterator it2 = it->begin(); it2 != it->end() ;it2++){
                int k = labels.at<int>(labels_i_w);
                
                if(k > 2){
                    k = 2;
                } 

                cell_class::immune mylabel = cell_class::setCellType( k );
                
                // assign center_value of region accordingly to the gmm-classificaton: single cells ; cluster ; other
                it2->setClass( mylabel );
                
                labels_i_w++;
                

                if (FLAG_DEBUG){
                    // get the center value of current ROI         
                    cv::Point center_point = it2->getCentroidPoint();
                    cell_class::immune myclass = it2->getClass();
                    std::cout << "\tcenter_point: " << center_point << "\tclass of ROI: " << cell_class::getCellClass(myclass) << "\t class of GMM: " << k << std::endl;
                }
                
            }

        }

    }

    /**
     * estimate Gaussian Mixture distribution of cell areas   
     */
    void estimate_gmm_cell_areas(const int &C, std::vector<std::vector<Region>> &regions, const int &DELAY, const std::string &INPUTDIR_COLOR, const std::string &OUTDIR, const bool &FLAG_DEBUG){

        cv::Ptr<cv::ml::EM> em_model_w = cv::ml::EM::create();        
        cv::Mat labels_w;

        Classification::approxMixDistrEM(C, regions, labels_w, em_model_w, OUTDIR, false);

        // save model
        std::string filename = OUTDIR + "/CvEM.yaml";
        em_model_w->save(filename);
        
        tracking_region_association::assign_gmm_class_to_region(regions, labels_w, FLAG_DEBUG); 
        
        if(FLAG_DEBUG){
			std::string file1 = OUTDIR + "/3_gmm_region_classification/";
            
            outputs::showClassifiedSegmentedRegions2D(file1, INPUTDIR_COLOR , DELAY, regions, true, false);
		}
        
    }

    /**
     * Assign only single cell class because in the data are no clusters.
     * This is because clusters are logical not possible in the data or removed during pre-processing.
     * 
     * @param regions all regions-per-image .
     * @param FLAG_DEBUG verbose parameter.
     */
    void assign_single_class_to_region(std::vector<std::vector<Region>> &regions, const bool &FLAG_DEBUG) {

        /// iterate over all images
        for(std::vector<std::vector<Region>>::iterator it = regions.begin(); it != regions.end(); it++){
            /// iterate over all regions
            for(std::vector<Region>::iterator it2 = it->begin(); it2 != it->end() ;it2++){
                
                // assign single cell class
                it2->setClass( cell_class::immune::SINGLE );
                
            }

        }
        
    }

} // tracking_region_association