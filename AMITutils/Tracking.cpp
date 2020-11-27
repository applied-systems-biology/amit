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

#include "Tracking.h"
#include "Calcs.h"
#include <lemon/list_graph.h>
#include <lemon/matching.h>
#include "RegionSplitting.h"


namespace Tracking
{

    /**
     *
     * @param current --> current region
     * @param regions --> vector of regions
     * @param images --> vector of images with regions ( of klass 1 for example)
     * @param t --> current time step
     * @param id --> current ID
     * @param track --> current track
     *
     * 1) go through all regions in the next time step (t+1)
     * 		- check klass of the region (has to be 1 or 2) (or just 1 for 2 classes)
     * 			- check if the region already has an ID and if it overlaps with the current region (at least 1 px)
     * 				- save region as successor
     * 				- increase number of overlaps
     *
     * 2) if ONE overlap was found and the successor has klass 1 and not already an ID:
     * 		- set ID of the successor region
     * 		- add successor region to track
     * 		- starting from here, find the next region to add
     *
     *
     */
    void findNext(Region &current, std::vector<std::vector<Region>> &regions, int t, int id, CellTrack &track){
        // counter for number of overlapping regions
        int n_overlaps = -1;

        /// map with < int = num_overlap_pixels, Region = corresponding overlap region >
        std::map<int, Region*> overlap_successors;

        // iterate through all regions of the next frame and look for overlaps
        for(auto & it : regions.at(t+1)){
            
            // check if one of the region pixels is white (use overlap image for faster computation)
            Region tempr = it;
                         
            // check for single cells and cell clusters
            if(tempr.getClass() == cell_class::immune::SINGLE || tempr.getClass() == cell_class::immune::CLUSTER){

                int number_of_overlapping_pixels = overlap(current, tempr);

                if (tempr.getId() == 0 && number_of_overlapping_pixels > 0) {
                    /// if overlap occur, store the regions with the number of overlapping pixels
                    overlap_successors[ number_of_overlapping_pixels ] = &it;

                    n_overlaps++;
                }

            }
        }

        /// extract region with the most overlapping regions in case there are more than one overlapping single cells
        Region *successor;
        if (overlap_successors.size() > 0) {
            int most_overlapping_regions = overlap_successors.rbegin()->first;
            successor = overlap_successors[most_overlapping_regions];
        }

        /// OLD: only exactly one overlap with a single cells
        /// NEW: if more than exactly one overlap with a single cells: assign sell with most overlapping pixels
        if(n_overlaps != -1 && successor->getClass() == cell_class::immune::SINGLE && successor->getId() == 0){

            // continue track if successor has class of single cell (0)
            successor->setId(id);

            // add to CellTrack
            track.add(*successor);

            // continue iteratively with next frame
            if((unsigned)t < regions.size() - 2){
                Tracking::findNext(*successor, regions, t+1, id, track);
            }
        }

    }

    /**
     *
     * @param current --> current region
     * @param regions --> vector of regions
     * @param images --> vector of images with regions ( of klass 1 for example)
     * @param t --> current time step
     * @param id --> current ID
     * @param track --> current track
     *
     * 1) go through all regions in the next time step (t+1)
     * 		- check klass of the region (has to be 1 or 2)
     * 			- check if the region already has an ID and if it overlaps with the current region (at least 1 px)
     * 				- save region as successor
     * 				- increase number of overlaps
     *
     * 2) if ONE overlap was found and the successor has klass 1 and not already an ID:
     * 		- set ID of the successor region
     * 		- add successor region to track
     * 		- starting from here, find the next region to add
     *
     *
     */
    void findNext(RegionP &current, std::vector<std::vector<RegionP>> &regions, std::vector<cv::Mat> &images, int t, int id, CellTrackP &track){
      
        /// counter for number of overlapping regions
        int n_overlaps = -1;

        /// map with < int = num_overlap_pixels, Region = corresponding overlap region >
        std::map<int, RegionP*> overlap_successors;

        // // iterate through all regions of the next frame and look for overlaps
        for(auto & it : regions.at(t+1)){
            // check if one of the region pixels is white (use overlap image for faster computation)
            RegionP tempr = it;

            /// deprecated: no other region types available        
            // check for single cells and cell clusters, CAUTION: assume: 0 -> single cell, 1 -> cell cluster, (2 -> other)
            
            int number_of_overlapping_pixels = overlap(current, tempr);

            if (tempr.getId() == 0 && number_of_overlapping_pixels > 0) {
                /// if overlap occur, store the regions with the number of overlapping pixels
                overlap_successors[ number_of_overlapping_pixels ] = &it;

                n_overlaps++;
            }

        }

        /// extract region with the most overlapping regions in case there are more than one overlapping single cells
        RegionP *successor;
        if (overlap_successors.size() > 0) {
            int most_overlapping_regions = overlap_successors.rbegin()->first;
            successor = overlap_successors[most_overlapping_regions];
        }

        /// OLD: only exactly one overlap with a single cells
        /// NEW: if more than exactly one overlap with a single cells: assign sell with most overlapping pixels
        if(n_overlaps != -1 && successor->getClass() == cell_class::immune::SINGLE && successor->getId() == 0){

            // continue track if successor has class of single cell (0)
            successor->setId(id);

            // add to CellTrack
            track.add(*successor);

            // continue iteratively with next frame
            if((unsigned)t < regions.size() - 2){
                Tracking::findNext(*successor, regions, images, t+1, id, track);
            }
        }

    }
    
    /**
     * This function performs nearest-neighbor association tracking using overlaps.
     * edited by PP to execute function with 2 and 3 classes
     *
     * @param regions input regions
     * @param images vector of images including only regions that can be tracked
     * @param tracks set of created cell tracks
     */
    void trackWithOverlap(std::vector<std::vector<Region>> &regions, std::vector<CellTrack> &tracks){

        std::vector<std::vector<Region>>::iterator vreg_it;
        std::vector<Region>::iterator reg_it;
        
        /// global ID-counter
        int id = 1; 

        int t = 0;

        /// determine the maximal number of frames
        int t_max = regions.size();

        /// start with t = 0; search in the following frames 
        for(vreg_it = regions.begin(); vreg_it != regions.end(); vreg_it++){
            for(reg_it = vreg_it->begin(); reg_it != vreg_it->end(); reg_it++){

                /// check for starting track and single cell - class
                if(reg_it->getId() == 0 && reg_it->getClass() == cell_class::immune::SINGLE) {
 
                    /// set ID
                    reg_it->setId(id);

                    /// create new CellTrack
                    CellTrack track(t, id);
                    track.add(*reg_it);

                    /// find successor regions in the following frames, if frame is NOT the last frame
                    if (t != t_max-1) {
                        Tracking::findNext(*reg_it, regions, t, id, track);
                    }

                    tracks.push_back(track);

                    id++;
                }
            }

            t++;
        }

    }

    /**
     *
     * @param regions --> pmn regions
     * @param images --> binary images of the regions (klass 1 only --> single cells)
     * @param tracks --> vector for resulting tracks
     *
     * 1) go through all frames:
     * 		- go through all regions in the frame:
     * 			- check, if the region doesn't already have an ID and is of klass 1 (single cell)
     * 				- set ID
     * 				- compute centroid
     * 				- start cell track --> create track-object and set first region
     * 				- find next region in the following frames
     * 					- only klass 1 regions with at least one pixel overlap
     * 					- only clear overlaps with exactly one other region are considered
     * 					- recursive approach
     * 				- save finished track
     *				- increment ID-counter
    *
    */
    void trackWithOverlap(std::vector<std::vector<RegionP>> &regions, std::vector<cv::Mat> images, std::vector<CellTrackP> &tracks){
        
        std::vector<std::vector<RegionP>>::iterator vreg_it;
        std::vector< RegionP>::iterator reg_it;
        int id = 1; //ID-Zähler
        int t = 0;

        /// start with t = 0; search in the following frames
        // for(vreg_it = regions.begin(); vreg_it != regions.end()-2; vreg_it++){
        for(vreg_it = regions.begin(); vreg_it != regions.end()-1; vreg_it++){ // TODO

            for(reg_it = vreg_it->begin(); reg_it != vreg_it->end(); reg_it++){

                if(reg_it->getId() == 0 && reg_it->getClass() == cell_class::immune::SINGLE ){
                    /// set ID
                    reg_it->setId(id);
                    reg_it->computeCentroid();

                    /// create new CellTrack
                    CellTrackP track(t, id);
                    track.add(*reg_it);

                    /// find successor regions in the following frames
                    findNext(*reg_it, regions, images, t, id, track);

                    tracks.push_back(track);

                    id++;
                }
            }

            t++;
        }

    }

    /**
     * This function creates Cell Tracks for Regions of class 1 (cell clusters).
     * All cell tracks have track length 1, so for each region one cell track is created.
     * edited by Philipp to execute function with 2 and 3 classes
     *
     * @param regions v<v<Region>> input regions
     */
    void trackWithOverlapClass2(std::vector<std::vector<Region>> &regions, std::vector<CellTrack> &tracks){
        
        std::vector<std::vector<Region>>::iterator vreg_it;
        std::vector<Region>::iterator reg_it;
       
        /// ID-counter
        int id = 10000; 
        int t = 0;

        /// start with t = 0; search in the following frames
        for(vreg_it = regions.begin(); vreg_it != regions.end(); vreg_it++){
            for(reg_it = vreg_it->begin(); reg_it != vreg_it->end(); reg_it++){
                
                /// check for starting track and cell cluster - class
                if( reg_it->getId() == 0 && reg_it->getClass() == cell_class::immune::CLUSTER ){
                
                    /// set ID
                    reg_it->setId(id);

                    /// create new CellTrack
                    CellTrack track(t, id);
                    track.add(*reg_it);

                    tracks.push_back(track);
                    id++;
                }
            }

            t++;
        }

    }

    /**
     * This function computes the standard deviation of cell areas for all cell tracks
     *
     * @param tracks1 set of cell tracks
     *
     * @return standard deviation of cell area
     */
    double getAreaDeviation(std::vector<CellTrack> &tracks1){
        double dev = 0.0;

        std::vector<double> areas;
        std::vector<CellTrack>::iterator trackit;
        
        for(trackit = tracks1.begin(); trackit != tracks1.end(); trackit++){
            areas.push_back(trackit->getAverageArea());
        }

        dev = Calcs::computeSd(areas);

        return dev;
    }

    /**
     * This function computes the standard deviation of cell areas for all cell tracks
     *
     * @param tracks1 set of cell tracks
     *
     * @return standard deviation of cell area
     */
    double getAreaDeviation(std::vector<CellTrackP> &tracks1){
        double dev = 0.0;

        std::vector<double> areas;
        
        std::vector<CellTrackP>::iterator trackit;
        for(trackit = tracks1.begin(); trackit != tracks1.end(); trackit++){
            areas.push_back((double)trackit->getAverageArea());
        }

        dev = Calcs::computeSd(areas);

        return dev;
    }

    /**
     *  @param tracks1 set of cell tracks
     *
     *  @return maximum speed
     */
    double getMaxSpeed(std::vector<CellTrack> &tracks1){
        double speed = 0.0;

        std::vector<CellTrack>::iterator trackit;
        
        for(trackit = tracks1.begin(); trackit != tracks1.end(); trackit++){
            double tmp = trackit->getMaxSpeed();
    
            if(tmp > speed){
                speed  = tmp;
            }
    
        }

        return speed;
    }

    /**
     *  @param tracks1 set of cell tracks
     *
     *  @return maximum speed
     */
    double getMaxSpeed(std::vector<CellTrackP> &tracks1){
        double speed = 0.0;

        std::vector<CellTrackP>::iterator trackit;
        
        for(trackit = tracks1.begin(); trackit != tracks1.end(); trackit++){
            double tmp = trackit->getMaxSpeed();
    
            if(tmp > speed){
                speed  = tmp;
            }
    
        }

        return speed;
    }

    /**
     * This function computes the total length of a set of cell tracks.
     *
     * @param tracks set of cell tracks
     *
     * @return total length
     */
    int getTracksLengthTotal(std::vector<CellTrack> &tracks){
        int length = 0;
        std::vector<CellTrack>::iterator trackit;

        for(trackit = tracks.begin(); trackit != tracks.end(); trackit++){
            length += trackit->getLength();
        }

        return length;
    }

    /**
     * This function computes the total length of a set of cell tracks.
     *
     * @param tracks set of cell tracks
     *
     * @return total length
     */
    int getTracksLengthTotal(std::vector<CellTrack*> &tracks){
        
        int length = 0;

        std::vector<CellTrack*>::iterator trackit;

        for(trackit = tracks.begin(); trackit != tracks.end(); trackit++){
            if(CellTrackP* p = dynamic_cast<CellTrackP*>(*trackit)){
                length += p->getLength();
            }
            else{
                length += (*trackit)->getLength();
            }
        }
        return length;
    }

    /**
     * This function adds a Region to a CellTrack given its ID at a specified time point.
     *
     * @param r Region
     * @param ct set of cell tracks
     * @param id cell track ID
     * @param t time point
     */
    void addRegionToTrack(Region &r, std::vector<CellTrack> &ct, const int &id, const int &t){

        std::vector<CellTrack>::iterator trackit;
        bool found =false;

        for(trackit = ct.begin(); trackit != ct.end(); trackit++){
            
            if(trackit->getID() == id){
                r.setId(id);
    
                // distinguish for (none) given timepoint
                 if (t != -1) {
                     trackit->add(r, t);
                 }
                 else {
                     trackit->add(r);
                 }

                found =  true;
                break;
            }
        }

        if(!found){
            std::cout << "Error: ID " << id << " not found (\"addRegionToTrack()\")" << std::endl;
        }
    }

    /**
     * This function adds a Region to a CellTrack given its ID at a specified time point.
     *
     * @param r Region
     * @param ct set of cell tracks
     * @param id cell track ID
     * @param t time point
     */
    void addRegionToTrack(RegionP &r, std::vector<CellTrackP> &ct, int id, int t){

        std::vector<CellTrackP>::iterator trackit;
        bool found =false;
        
        for(trackit = ct.begin(); trackit != ct.end(); trackit++){
            if(trackit->getID() == id){
                r.setId(id);
                trackit->add(r, t);
                found =  true;
                break;
            }
        }
        if(!found){
            std::cout << "Error: ID " << id << " not found (\"addRegionToTrack()\")" << std::endl;
        }

    }

    /**
     * This function adds the IDs of interacting cells to a CellTrack given its ID
     *
     * @param tracks set of cell track in which trackid exists
     * @param trackid ID of cell track
     * @param t time point of interaction
     * @interacting_cell set of interacting cells
     */
    void addInteractionToTrack(std::vector<CellTrack> &tracks, const int trackid, const int t, std::vector<Region> &interacting_cells){

        bool found = false;

        std::vector<CellTrack>::iterator trackit;

        for(trackit = tracks.begin(); trackit != tracks.end(); trackit++){
            if(trackit->getID() == trackid){

                for(int i = 0; i < interacting_cells.size(); i++){
                    if(interacting_cells.at(i).getId() != trackid){
                        trackit->setInteraction(t, interacting_cells.at(i).getId());
                    }
                }

                found = true;
                break;
            }
        }

        if(!found){
            std::cout << "Error: ID " << trackid << " not found (\"addInteractionToTrack()\")" << std::endl;
        }
    }

    /**
     * This function adds the IDs of interacting cells to a CellTrack givent its ID
     *
     * @param tracks set of cell track in which trackid exists
     * @param trackid ID of cell track
     * @param t time point of interaction
     * @interacting_cell set of interacting cells
     */
    void addInteractionToTrack(std::vector<CellTrackP> &tracks, const int trackid, const int t, std::vector<RegionP> &interacting_cells){

        bool found =false;

        std::vector<CellTrackP>::iterator trackit;

        for(trackit = tracks.begin(); trackit != tracks.end(); trackit++){
            if(trackit->getID() == trackid){
    
                for(size_t i = 0; i < interacting_cells.size(); i++){
                    if(interacting_cells.at(i).getId() != trackid){
                        trackit->setInteraction(t, interacting_cells.at(i).getId());
                    }
                }

                found = true;
                break;
            }
        }

        if(!found){
            std::cout << "Error: ID " << trackid << " not found (\"addInteractionToTrack()\")" << std::endl;
        }
    }

    /**
     * This function deletes a region specified by its ID in a two-dimensional set of regions at the given time point.
     *
     * @param regions two-dimensional set of regions
     * @param t time point
     * @param ID region ID
     */
    void deleteRegion(std::vector<std::vector<Region>> &regions, const int &t, const int &ID){
        std::vector<Region>::iterator regit;

        regit = regions.at(t).begin();
        std::vector<Region> tmp;

        bool found = false;
        for(; regit != regions.at(t).end(); regit++){
            if(regit->getId() != ID){
                tmp.push_back(*regit);
            }
            else{
                found = true;
            }
        }

        regions.at(t) = tmp;

        if(!found){
            std::cout << "ERROR: DELETE REGION, ID"  << ID << "( " << t << ") NOT FOUND!!!" << std::endl;
        }

    }

    /**
     * This function deletes a region specified by its ID in a two-dimensional set of regions at the given time point.
     *
     * @param regions two-dimensional set of regions
     * @param t time point
     * @param ID region ID
     */
    void deleteRegion(std::vector<std::vector<RegionP>> &regions, int t, int ID){
        
        std::vector<RegionP>::iterator regit;
        regit = regions.at(t).begin();
        
        std::vector<RegionP> tmp;

        bool found = false;
        for(; regit != regions.at(t).end(); regit++){
            if(regit->getId() != ID){
                tmp.push_back(*regit);
            }
            else{
                found = true;
            }
        }

        regions.at(t) = tmp;

        if(!found){
            std::cout << "ERROR: DELETE REGION, ID"  << ID << "( " << t << ") NOT FOUND!!!" << std::endl;
        }

    }

    /**
     *	This function reduces the set of overlapping regions by its smallest element.
    *
    *	@param cluster clustered region
    *	@param overlaps set of overlapping regions
    *
    */
    void reduceOverlaps(Region &cluster, std::vector<Region> &overlaps){
        std::vector<Region>::iterator regit;
        std::vector<int> overlappercell;
        for(regit =  overlaps.begin(); regit != overlaps.end(); regit++){
            int o = overlapN(cluster, *regit);
            overlappercell.push_back(o);
        }

        // remove smallest overlapregion
        int i = 0;
        int min = INT_MAX;
        std::vector<int>::iterator intit;
        for(size_t j = 0; j < overlappercell.size(); j++){
            if(overlappercell.at(j) < min){
                i = j;
                min = overlappercell.at(j);
            }
        }

        // remove ith (smallest) element of vector overlaps
        std::vector<Region> tmp;
        for(size_t j = 0; j < overlaps.size(); j++){
            if((signed) j != i){
                tmp.push_back(overlaps.at(j));
            }
        }

        overlaps = tmp;
    }

    /**
     *	This function reduces the set of overlapping regions by its smallest element.
    *
    *	@param cluster clustered region
    *	@param overlaps set of overlapping regions
    *
    */
    void reduceOverlaps(Region *cluster, std::vector<Region*> &overlaps){
        
        std::vector<Region*>::iterator regit;
        std::vector<int> overlappercell;
        
        for(regit =  overlaps.begin(); regit != overlaps.end(); regit++){
            int o = overlapN(*cluster, **regit);
            overlappercell.push_back(o);
        }

        /// remove smallest overlapregion
        int i = 0;
        int min = INT_MAX;
        
        std::vector<int>::iterator intit;
        for(size_t j = 0; j < overlappercell.size(); j++){
            if(overlappercell.at(j) < min){
                i = j;
                min = overlappercell.at(j);
            }
        }

        /// remove ith (smallest) element of vector overlaps
        std::vector<Region*> tmp;
        for(size_t  j = 0; j < overlaps.size(); j++){
            if((signed) j != i){
                tmp.push_back(overlaps.at(j));
            }
        }
        
        overlaps = tmp;

    }

    void computeOverlapRegionsForCellTracks(std::vector<CellTrack> &tracks, std::vector<CellTrack> &tracks_class2, int t_max){
        std::vector< CellTrack>::iterator track_it, track2_it;

        for(track2_it =  tracks_class2.begin(); track2_it != tracks_class2.end(); track2_it++){
            track2_it->clearOverlaps();

            int t_end = track2_it->getEndTime();
            int t_start = track2_it->getStartTime();

            Region tmp_start, tmp_end;
            track2_it->getFirstRegion(tmp_start);
            track2_it->getLastRegion(tmp_end);

            //look forward
            //compute overlap regions for class 1
            if(t_end < t_max-1){
                for(track_it = tracks.begin(); track_it != tracks.end(); track_it++){
                    if(track_it->getStartTime() == t_end+1){
                        Region tmp;
                        track_it->getRegionAt(t_end+1, tmp);
                        if(overlap(tmp_end, tmp)){
                            track2_it->overlapregionsf1.push_back(tmp);
                        }
                    }
                }

                //compute overlap regions for class 2
                for(track_it = tracks_class2.begin(); track_it != tracks_class2.end(); track_it++){
                    if(track_it->getStartTime() == t_end+1){
                        Region tmp;
                        track_it->getRegionAt(t_end+1, tmp);
                        if(overlap(tmp_end, tmp)){
                            track2_it->overlapregionsf2.push_back(tmp);
                        }
                    }
                }
            }


            //look backward
            //compute overlap regions for class 1
            if(t_start > 0){
                for(track_it = tracks.begin(); track_it != tracks.end(); track_it++){
                    if(track_it->getEndTime() == t_start -1){
                        Region tmp;
                        track_it->getRegionAt(t_start-1, tmp);
                        if(overlap(tmp_start, tmp)){
                            track2_it->overlapregionsb1.push_back(tmp);
                        }
                    }
                }

                //compute overlap regions for class 2
                for(track_it = tracks_class2.begin(); track_it != tracks_class2.end(); track_it++){
                    if(track_it->getEndTime() == t_start -1){
                        Region tmp;
                        track_it->getRegionAt(t_start-1, tmp);
                        if(overlap(tmp_start, tmp)){
                            track2_it->overlapregionsb2.push_back(tmp);
                        }
                    }
                }

            }


        }
    }

    void computeOverlapRegionsForCellTracks(std::vector<CellTrackP> &tracks, std::vector<CellTrackP> &tracks_class2, int t_max){
        
        std::vector<CellTrackP>::iterator track_it, track2_it;

        for(track2_it =  tracks_class2.begin(); track2_it != tracks_class2.end(); track2_it++){
    
            track2_it->clearOverlaps();

            std::vector<RegionP> track;
            track2_it->getTrack(track);

            int t_end = track2_it->getEndTime();
            int t_start = track2_it->getStartTime();

            RegionP tmp_start, tmp_end;
            track2_it->getFirstRegion(tmp_start);
            track2_it->getLastRegion(tmp_end);

            /// look forward
            // compute overlap regions for class 1
            if(t_end < t_max-1){
                
                for(track_it = tracks.begin(); track_it != tracks.end(); track_it++){
                    if(track_it->getStartTime() == t_end+1){
                        RegionP tmp;
                        track_it->getRegionAt(t_end+1, tmp);
                
                        if(overlap(tmp_end, tmp)){
                            track2_it->overlapregionsf1.push_back(tmp);
                        }
                    }
                }

                // compute overlap regions for class 2
                for(track_it = tracks_class2.begin(); track_it != tracks_class2.end(); track_it++){
                    if(track_it->getStartTime() == t_end+1){
                        RegionP tmp;
                        track_it->getRegionAt(t_end+1, tmp);
    
                        if(overlap(tmp_end, tmp)){
                            track2_it->overlapregionsf2.push_back(tmp);
                        }
                    }
                }
            }

            /// look backward
            // compute overlap regions for class 1
            if(t_start > 0){
               
                for(track_it = tracks.begin(); track_it != tracks.end(); track_it++){
                    if(track_it->getEndTime() == t_start -1){
                        RegionP tmp;
                        track_it->getRegionAt(t_start-1, tmp);
               
                        if(overlap(tmp_start, tmp)){
                            track2_it->overlapregionsb1.push_back(tmp);
                        }
                    }
                }

                // compute overlap regions for class 2
                for(track_it = tracks_class2.begin(); track_it != tracks_class2.end(); track_it++){
                    if(track_it->getEndTime() == t_start -1){
                        RegionP tmp;
                        track_it->getRegionAt(t_start-1, tmp);
    
                        if(overlap(tmp_start, tmp)){
                            track2_it->overlapregionsb2.push_back(tmp);
                        }
                    }
                }
            }

        }
    }

    void computeOverlapRegionsForCellTracksFC(std::vector<CellTrackP> &tracks_fc, std::vector<CellTrackP> &tracks_class2_fc, std::vector<CellTrackP> &tracks_class2, int t_max){
        
        std::vector< CellTrackP>::iterator track_it, track2_it;

        for(track2_it =  tracks_class2.begin(); track2_it != tracks_class2.end(); track2_it++){
            
            track2_it->clearOverlaps();

            std::vector<RegionP> track;
            track2_it->getTrack(track);

            int t_end = track2_it->getEndTime();
            int t_start = track2_it->getStartTime();

            RegionP tmp_start, tmp_end;
            track2_it->getFirstRegion(tmp_start);
            track2_it->getLastRegion(tmp_end);

            /// look forward
            // compute overlap regions for class 1 FC tracks
            if(t_end < t_max-1){
                for(track_it = tracks_fc.begin(); track_it != tracks_fc.end(); track_it++){
                    if(track_it->getStartTime() == t_end+1){
                        RegionP tmp;
                        track_it->getRegionAt(t_end+1, tmp);
                        if(overlap(tmp_end, tmp)){
                            track2_it->overlapregionsf1_fc.push_back(tmp);
                        }
                    }
                }

                // compute overlap regions for class 2
                for(track_it = tracks_class2_fc.begin(); track_it != tracks_class2_fc.end(); track_it++){
                    if(track_it->getStartTime() == t_end+1){
                        RegionP tmp;
                        track_it->getRegionAt(t_end+1, tmp);
                        if(overlap(tmp_end, tmp)){
                            track2_it->overlapregionsf2_fc.push_back(tmp);
                        }
                    }
                }
            }

            /// look backward
            // compute overlap regions for class 1
            if(t_start > 0){
                for(track_it = tracks_fc.begin(); track_it != tracks_fc.end(); track_it++){
                    if(track_it->getEndTime() == t_start -1){
                        RegionP tmp;
                        track_it->getRegionAt(t_start-1, tmp);
               
                        if(overlap(tmp_start, tmp)){
                            track2_it->overlapregionsb1_fc.push_back(tmp);
                        }
                    }
                }

                /// compute overlap regions for class 2
                for(track_it = tracks_class2_fc.begin(); track_it != tracks_class2_fc.end(); track_it++){
                    if(track_it->getEndTime() == t_start -1){
                        RegionP tmp;
                        track_it->getRegionAt(t_start-1, tmp);
                        if(overlap(tmp_start, tmp)){
                            track2_it->overlapregionsb2_fc.push_back(tmp);
                        }
                    }
                }
            }
        }
    }

    void computeOverlapRegionsForCellTrack(CellTrackP &track, std::vector<CellTrackP> &tracks, std::vector<CellTrackP> &tracks_class2, int t_max){
	    
        std::vector<CellTrackP>::iterator track_it;

		track.clearOverlaps();

		int t_end = track.getEndTime();
		int t_start = track.getStartTime();

		RegionP tmp_start, tmp_end;
		track.getFirstRegion(tmp_start);
		track.getLastRegion(tmp_end);

		/// look forward 
        // compute overlap regions for class 1
		if(t_end < t_max-1){
			for(track_it = tracks.begin(); track_it != tracks.end(); track_it++){
				if(track_it->getStartTime() == t_end+1){
					RegionP tmp;
					track_it->getRegionAt(t_end+1, tmp);
					if(overlap(tmp_end, tmp)){
						track.overlapregionsf1.push_back(tmp);
                    }
				}
			}

			// compute overlap regions for class 2
			for(track_it = tracks_class2.begin(); track_it != tracks_class2.end(); track_it++){
				if(track_it->getStartTime() == t_end+1){
					RegionP tmp;
					track_it->getRegionAt(t_end+1, tmp);
					if(overlap(tmp_end, tmp)){
						track.overlapregionsf2.push_back(tmp);
					}
				}
			}
		}

		/// look backward
		// compute overlap regions for class 1
		if(t_start > 0){
			for(track_it = tracks.begin(); track_it != tracks.end(); track_it++){
				if(track_it->getEndTime() == t_start -1){
					RegionP tmp;
					track_it->getRegionAt(t_start-1, tmp);
			
            		if(overlap(tmp_start, tmp)){
						track.overlapregionsb1.push_back(tmp);
					}
				}
			}

			// compute overlap regions for class 2
			for(track_it = tracks_class2.begin(); track_it != tracks_class2.end(); track_it++){
				if(track_it->getEndTime() == t_start -1){
					RegionP tmp;
					track_it->getRegionAt(t_start-1, tmp);
	
    				if(overlap(tmp_start, tmp)){
						track.overlapregionsb2.push_back(tmp);
					}
				}
			}
		}

    }

    /**
     * This function splits a cluster dependent on its area.
     * The number of cells in the cluster are estimated by the typical area of a cell. 
     * If the area criterion (|sum(A_overlaps) - A_cluster| < SD_area * c * n_overlaps) does not hold, the size of overlapping regions is subsequently reduced.
     *
     * @param t time point
     * @param areadev standard deviation of cell areas
     * @param c
     * @param overlap set of overlapping regions
     * @param tracks set of cell tracks of class 1
     * @param track2 current class-2-cell-track
     * @param regions two-dimensional set of regions
     */
    bool splitClusterDependentOnArea(int t, double areadev, double c, std::vector<Region> &overlap, std::vector<CellTrack> &tracks, CellTrack &track2, std::vector<std::vector<Region>> &regions){

        int area = 0;
        
        std::vector<Region>::iterator regit;
        for(regit = overlap.begin(); regit != overlap.end(); regit++){
            area += regit->getArea();
        }

        Region cluster;
        track2.getRegionAt(t, cluster);

        if( abs(area - cluster.getArea()) < areadev * c * overlap.size() ){
                RegionSplitting::separateClusterGMM(cluster, overlap);

                // add regions to cell tracks
                for(size_t i = 0; i < overlap.size(); i++){
                    Tracking::addRegionToTrack(overlap.at(i), tracks, overlap.at(i).getId(), t);

                    overlap.at(i).setClass( cell_class::immune::SINGLE );
                    
                    regions.at(t).push_back(overlap.at(i));
                    if(overlap.size() > 1){
                        Tracking::addInteractionToTrack(tracks, overlap.at(i).getId(), t, overlap);
                    }
                }

                // remove cluster region from regions vector
                Tracking::deleteRegion(regions, t, track2.getID());

                return true;
        }
        // area criterion does not hold -> check the overlapping regions
        else{

            // cluster consists of at least two regions AND the sum of the subregions is bigger than the cluster
            if(overlap.size() > 1 && area > cluster.getArea()){ 
                reduceOverlaps(cluster, overlap);

                if(Tracking::splitClusterDependentOnArea(t, areadev, c, overlap, tracks, track2, regions)){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        }

    }

    /**
     * This function includes calls 2 regions into the set of cell tracks by splitting the clusters and associating the subregions to cell tracks.
     *
     * @param regions two-dimensional set of regions
     * @param tracks set of cell tracks class 1
     * @param tracks_class2 set of cell tracks class 2
     * @param areadev standard deviation of cell areas
     * @param c
     */
    void trackInteractingRegions(std::vector<std::vector<Region>> &regions, std::vector<CellTrack> &tracks, std::vector<CellTrack> &tracks_class2, double areadev, double c){
        std::vector< CellTrack>::iterator track_it, track2_it;
        int t_max = regions.size();

        std::vector<CellTrack> track_tmp2;

        Tracking::computeOverlapRegionsForCellTracks(tracks, tracks_class2, t_max);

        for(track2_it =  tracks_class2.begin(); track2_it != tracks_class2.end(); track2_it++){

            int t = track2_it->getStartTime();

            std::vector<Region> overlapb, overlapf;

            int b1 = track2_it->overlapregionsb1.size();
            int b2 = track2_it->overlapregionsb2.size();
            int f1 = track2_it->overlapregionsf1.size();
            int f2 = track2_it->overlapregionsf2.size();

            //create vectors for complete bw overlaps and fw overlaps
            overlapb = track2_it->overlapregionsb1;
            overlapb.insert(overlapb.end(), track2_it->overlapregionsb2.begin(), track2_it->overlapregionsb2.end());
            overlapf = track2_it->overlapregionsf1;
            overlapf.insert(overlapf.end(), track2_it->overlapregionsf2.begin(), track2_it->overlapregionsf2.end());


            if(b1 == 0 && b2 == 0 && f1 == 0 && f2 == 0){
                // remove track
                // continue;
            }

            else if(b1 == 0 && b2 == 0 && f1 > 0 && f2 == 0){ //0,0,n>1,0
                std::vector<Region> overlap = track2_it->overlapregionsf1;

                Region cluster;
                track2_it->getRegionAt(t, cluster);

                RegionSplitting::separateClusterGMM(cluster, overlap);

                // add regions to cell tracks
                for(int i = 0; i < overlap.size(); i++){
                    Tracking::addRegionToTrack(overlap.at(i), tracks, overlap.at(i).getId(), t);

                    overlap.at(i).setClass( cell_class::immune::SINGLE );

                    regions.at(t).push_back(overlap.at(i));
                    if(overlap.size() > 1){
                        Tracking::addInteractionToTrack(tracks, overlap.at(i).getId(), t, overlap);
                    }
                }

                //remove cluster region from regions vector
                Tracking::deleteRegion(regions, t, track2_it->getID());

            }

            else if(b1 > 0 && b2 == 0 && f1 == 0 && f2 == 0){ //n>0,0,0,0
                std::vector<Region> overlap = track2_it->overlapregionsb1;

                Region cluster;
                track2_it->getRegionAt(t, cluster);

                RegionSplitting::separateClusterGMM(cluster, overlap);

                //add regions to cell tracks
                for(int i = 0; i < overlap.size(); i++){
                    Tracking::addRegionToTrack(overlap.at(i), tracks, overlap.at(i).getId());

                    overlap.at(i).setClass( cell_class::immune::SINGLE );

                    regions.at(t).push_back(overlap.at(i));
                    if(overlap.size() > 1){
                        Tracking::addInteractionToTrack(tracks, overlap.at(i).getId(), t, overlap);
                    }
                }

                //remove cluster region from regions vector
                Tracking::deleteRegion(regions, t, track2_it->getID());
            }

            else if(b1 > 0 && b1 == f1 && b2 == 0 && f2 == 0){ //n,0,n,0
                //split region
                std::vector<Region> overlap = track2_it->overlapregionsb1;

                Region cluster;
                track2_it->getRegionAt(t, cluster);

                RegionSplitting::separateClusterGMM(cluster, overlap);

                //add regions to cell tracks
                for(int i = 0; i < overlap.size(); i++){
                    Tracking::addRegionToTrack(overlap.at(i), tracks, overlap.at(i).getId());

                    overlap.at(i).setClass( cell_class::immune::SINGLE ); 

                    regions.at(t).push_back(overlap.at(i));
                    if(overlap.size() > 1){
                        Tracking::addInteractionToTrack(tracks, overlap.at(i).getId(), t, overlap);
                    }
                }

                //remove cluster region from regions vector
                Tracking::deleteRegion(regions, t, track2_it->getID());

            }

            else if(b1 > 1 && b2 == 0 && f1 >= 0 && f2 >= 1){ //n>1,0,m>=0,1
                //split cluster into first n regions
                std::vector<Region> overlap = track2_it->overlapregionsb1;

                if( ! Tracking::splitClusterDependentOnArea(t, areadev, c, overlap, tracks, *track2_it, regions)){
                    track_tmp2.push_back(*track2_it);
                }
                else{
//				cout << "splitted" << endl;
                }

            }


            else if(b1 >= 0 && b2 >= 1 && f1 > 1 && f2 == 0){ //m>=0,1,n>1,0
                //split cluster into last n regions
                std::vector<Region> overlap = track2_it->overlapregionsf1;

                if( ! Tracking::splitClusterDependentOnArea(t, areadev, c, overlap, tracks, *track2_it, regions)){
                    track_tmp2.push_back(*track2_it);
                }
                else{
//				cout << "splitted" << endl;
                }

            }

            else if(b1 == 0 && b2 >= 1 && f1 >= 1 && f2 == 0){ //0,>=1,>=1,0
                //split cluster into next n regions if area criterium holds
                std::vector<Region> overlap = track2_it->overlapregionsf1;

                if( ! Tracking::splitClusterDependentOnArea(t, areadev, c, overlap, tracks, *track2_it, regions)){
                    track_tmp2.push_back(*track2_it);
                }
                else{
//				cout << "splitted" << endl;
                }
            }

            else if(b1 >= 1 && b2 == 0 && f1 == 0 && f2 >= 1){ //>=1,0,0,1
                //split cluster into next n regions if area criterium holds
                std::vector<Region> overlap = track2_it->overlapregionsb1;


                if( ! Tracking::splitClusterDependentOnArea(t, areadev, c, overlap, tracks, *track2_it, regions)){
                    track_tmp2.push_back(*track2_it);
                }
                else{
//				cout << "splitted" << endl;
                }

            }

                //keiner der Fälle trifft zu
            else { //keep class 2 track
                track_tmp2.push_back(*track2_it);
            }

        }

        // ersetze class2-tracks durch vector ohne gesplittete tracks
        tracks_class2 = track_tmp2;

    }

    void createKlassImages(std::vector<cv::Mat> &images_klass0, std::vector<cv::Mat> &images_klass1, std::vector<cv::Mat> &images_klass2, std::vector<std::vector<RegionP>> &regions_over_time, int cols, int rows){
        
        images_klass0.clear();
        images_klass1.clear();
        images_klass2.clear();

        std::vector<std::vector<RegionP>>::iterator vreg_it;
        std::vector<RegionP>::iterator reg_it;
        std::vector<cv::Mat>::iterator it0, it1, it2;

        /// fill images with regions of specific class
        for(vreg_it = regions_over_time.begin(); vreg_it != regions_over_time.end(); vreg_it++){
            
            cv::Mat im0(rows, cols, CV_8UC1, cv::Scalar(0));
            cv::Mat im1(rows, cols, CV_8UC1, cv::Scalar(0));
            cv::Mat im2(rows, cols, CV_8UC1, cv::Scalar(0));
            
            for(reg_it = vreg_it->begin(); reg_it != vreg_it->end(); reg_it++){

                cell_class::immune k = reg_it->getClass();
                
                if (k == cell_class::immune::CLUSTER ) {
                    drawRegionOnImage(im2, (*reg_it));
                }
                else if(k == cell_class::immune::SINGLE ){ //single cell
                    drawRegionOnImage(im1, (*reg_it));
                }
                else{ //noise
                    drawRegionOnImage(im0, (*reg_it));
                }
            }

    //		cv::imshow("im1", im1);
    //		cv::waitKey(0);

            images_klass0.push_back(im0);
            images_klass1.push_back(im1);
            images_klass2.push_back(im2);
        }
    }

    /**
     * This function draws the region on the image (8bit image, white region)
     *
     * @param image input image
     * @param region region
     */
    void drawRegionOnImage(cv::Mat &image, const Region &region){
        
        std::vector<cv::Vec2i>::const_iterator pixit;

        if(image.channels() == 1){
                for(pixit = region.region_pixels.begin(); pixit != region.region_pixels.end(); pixit++){
                cv::Point p(pixit->val[1], pixit->val[0]);
                image.at<uchar>(p) = 255;
            }
        }
        else{
            std::cout << "Error: drawRegionOnImage(), channel != 1" << std::endl;
        }
    }

    void combineCellTracksWithGapSize(std::vector<CellTrack> &tracks1, const int t_max, const int gapsize, const int maxdistance, int rows, int cols, const bool &DEBUG){

        std::vector<CellTrack>::iterator trackit;
        std::vector<cv::Point> points;

        for(int t = t_max; t >= 0 ; t--){
            std::vector< int > id_start, id_start2;
            std::vector< int > id_end, id_end2;
            std::vector< int > gapsizes;
            std::vector<Region> reg_start, reg_end, reg_start2, reg_end2;

            for(trackit =  tracks1.begin(); trackit != tracks1.end(); trackit++){

                int t_end = trackit->getEndTime();
                int t_start = trackit->getStartTime();

                if(t_end == t){
                    id_end.push_back(trackit->getID());
                    Region tmp;
    
                    trackit->getLastRegion(tmp);
                    reg_end.push_back(tmp);

                    points.push_back(tmp.getCentroidPoint());
                }

                if(gapsize >= 0){

                    if(t_start == t+gapsize+1){

                        id_start.push_back(trackit->getID());
                        Region tmp;

                        trackit->getFirstRegion(tmp);
                        reg_start.push_back(tmp);

                        int d = t_start-t-1;
                        gapsizes.push_back(d);

                        points.push_back(tmp.getCentroidPoint());
                    }
                }
                /// negative gapsize
                else{ 
                    
                    if(t_start == t+gapsize+1){
                        id_start.push_back(trackit->getID());
                        Region tmp;
                    
                        trackit->getFirstRegion(tmp);
                        reg_start.push_back(tmp);

                        int d = t_start-t-1;
                        gapsizes.push_back(d);

                        points.push_back(tmp.getCentroidPoint());
                    }
                }
            }

            if(id_start.size() > 0 && id_end.size() > 0){
                combineTracksWithGraph(tracks1, id_start, id_end, reg_start, reg_end,  gapsizes, maxdistance, rows, cols, DEBUG);
            }

        }

    }

    /**
     *
     */
    void combineTracksWithGraph(std::vector<CellTrack> &tracks1, std::vector<int> &id_start, std::vector<int> &id_end, std::vector<Region> &reg_start, std::vector<Region> &reg_end, std::vector<int> gapsizes, double maxspeed, int rows, int cols, const bool &DEBUG){

        cv::Mat dist(id_start.size(), id_end.size(), CV_32F, cv::Scalar(0));
        Tracking::computeDistances(dist, id_start, id_end, reg_start, reg_end);

        cv::Mat dist2(id_start.size(), id_end.size(), CV_32F, cv::Scalar(0));
        Tracking::computeDistancesToImageBoundaries(dist2, id_start, id_end, reg_start, reg_end, rows, cols);
    
        /// solve with graph theory and declare graph object
        lemon::ListGraph g; 

        /// add nodes to the graph
        int n_nodes = id_end.size() + id_start.size();
        int n_edges = id_end.size() * id_start.size();

        for(int i = 0; i <  n_nodes; i++){
            g.addNode();
        }

        for(size_t i = 0; i < id_end.size(); i++){
            for(int j = id_end.size(); j < n_nodes; j++){
                g.addEdge(g.nodeFromId(i), g.nodeFromId(j));
            }
        }

        double min, max;
        cv::minMaxLoc(dist, &min, &max);
    
        /// add weights to an EdgeMap (type: float)
        lemon::ListGraph::EdgeMap<float> cost(g);
        
        int index = 0;
        for(int i = 0; i < dist.cols ; i++){
            for(int j = 0; j < dist.rows ; j++){
    
                int gap = std::abs(gapsizes.at(j));
    
                if( dist.at<float>(j,i) > (gap+1)*maxspeed){
                    cost[g.edgeFromId(index)] = minFloat;
                }
                else if(dist.at<float>(j,i) > dist2.at<float>(j,i)){
                    cost[g.edgeFromId(index)] = minFloat;
                }
                else{
                    cost[g.edgeFromId(index)] = (float) max + 1.f - dist.at<float>(j,i); //take negative weights to determine the minimum weight matching
                }
    
                index++;
            }
        }

        /// create object for maximum weighted graph matching
        lemon::MaxWeightedMatching<lemon::ListGraph, lemon::ListGraph::EdgeMap<float> > mwm(g, cost);
        /// run the algorithm
        mwm.run(); 

        std::vector<int> match;
        /// return the matchings

        int id1 = 0, id2 = 0;
        for(int i = 0; i < n_edges; i++){
            
            lemon::ListGraph::Edge e = g.edgeFromId(i);

            if (DEBUG) {
                std::cout << i << ": " << mwm.matching(e) << std::endl;
            }

            if(mwm.matching(e)){
                if(cost[e] > 0){
                    
                    if (DEBUG) {
                        std::cout << "match " << id_end.at(id1) << ", " << id_start.at(id2) << std::endl;
                    }
                    
                    Tracking::combineTwoCellTracks(tracks1, id_end.at(id1), id_start.at(id2));
                }
            }

            if (DEBUG) {
                std::cout << "id1 " << id_end.at(id1) << ", id2 " <<  id_start.at(id2) << std::endl;
            }
            
            id2++;

            /// increase id1-counter in every xth step (x=#tracks with starting time within the gapsize)
            if(i > 0 && ((i+1) % id_start.size()) == 0){
                id1++;
            }
            else if(i == 0 && id_start.size( ) == 1){
                id1++;
            }

            /// reset id2-counter in every xth step (x=#tracks with starting time within the gapsize)
            if(id2 % id_start.size() == 0){
                id2 = 0;
            }
        }
    }

    void combineTracksWithGraph(std::vector<CellTrackP> &tracks1, std::vector<int> & id_start, std::vector<int> &id_end, std::vector<RegionP> &reg_start, std::vector<RegionP> &reg_end, std::vector<int> gapsizes, double maxspeed){

        cv::Mat dist(id_start.size(), id_end.size(), CV_32F, cv::Scalar(0));
        Tracking::computeDistances(dist, id_start, id_end, reg_start, reg_end);

        /// solve with graph theory and declare graph object
        lemon::ListGraph g; 

        //// add nodes to the graph
        int n_nodes = id_end.size() + id_start.size();
        int n_edges = id_end.size() * id_start.size();

        for(int i = 0; i <  n_nodes; i++){
            g.addNode();
        }

        for(size_t i = 0; i < id_end.size(); i++){
            for(int j = id_end.size(); j < n_nodes; j++){
                g.addEdge(g.nodeFromId(i), g.nodeFromId(j));
            }
        }

        double min, max;
        cv::minMaxLoc(dist, &min, &max);
    
        /// add weights to an EdgeMap (type: float)
        lemon::ListGraph::EdgeMap<float> cost(g);
        int index = 0;
        
        for(int i = 0; i < dist.cols ; i++){
            for(int j = 0; j < dist.rows ; j++){
                int gap = std::abs(gapsizes.at(j));
    
                if( dist.at<float>(j,i) > (gap+1)*maxspeed){
                    cost[g.edgeFromId(index)] = minFloat;
                }
                else{
                    cost[g.edgeFromId(index)] = (float) max + 1.f - dist.at<float>(j,i); //take negative weights to determine the minimum weight matching
                }
    
                index++;
            }
        }

        /// create object for maximum weighted graph matching
        lemon::MaxWeightedMatching<lemon::ListGraph, lemon::ListGraph::EdgeMap<float> > mwm(g, cost);
        /// run the algorithm
        mwm.run(); 

        std::vector<int> match;
        /// return the matchings

        int id1 = 0, id2 = 0;
        for(int i = 0; i < n_edges; i++){
            lemon::ListGraph::Edge e = g.edgeFromId(i);

            if(mwm.matching(e)){
                if(cost[e] > 0){
                    Tracking::combineTwoCellTracks(tracks1, id_end.at(id1), id_start.at(id2));
                }
            }

            id2++;

            /// increase id1-counter in every xth step (x=#tracks with starting time within the gapsize)
            if(i > 0 && ((i+1) % id_start.size()) == 0){
                id1++;
            }
            else if(i == 0 && id_start.size( ) == 1){
                id1++;
            }

            /// reset id2-counter in every xth step (x=#tracks with starting time within the gapsize)
            if(id2 % id_start.size() == 0){
                id2 = 0;
            }
        }
    }

    void combineTracksWithGraph2(std::vector<CellTrackP> &tracks1, std::vector<int> & id_start, std::vector<int> & id_end, std::vector<RegionP> & reg_start, std::vector<RegionP> & reg_end, std::vector<int> gapsizes, double maxspeed){

        std::string s;

        /// at least one cell that has phagocytozed --> compute distance with Phagocytosis
        bool withP = false;

        for(size_t i = 0; i < id_start.size(); i++){
            for(size_t j = 0; j < id_end.size(); j++){
                if(reg_start.at(i).getIState() >= 3 || reg_end.at(j).getIState() >= 3){
                    withP = true;
                }
            }
        }

        cv::Mat dist(id_start.size(), id_end.size(), CV_32F, cv::Scalar(0));
        if(withP){
            Tracking::generateCostFunctionPairwiseEuklideanDistanceAndPhagocytosis(dist, reg_start, reg_end);
        }
        else{
            Tracking::computeDistances(dist, id_start, id_end, reg_start, reg_end);
        }

        /// solve with graph theory
        lemon::ListGraph g; //declare graph object

        /// add nodes to the graph
        int n_nodes = id_end.size() + id_start.size();
        int n_edges = id_end.size() * id_start.size();

        for(int i = 0; i <  n_nodes; i++){
            g.addNode();
        }

        for(size_t i = 0; i < id_end.size(); i++){
            for(int j = id_end.size(); j < n_nodes; j++){
                g.addEdge(g.nodeFromId(i), g.nodeFromId(j));
            }
        }

        double min, max;
        cv::minMaxLoc(dist, &min, &max);
    
        /// add weights to an EdgeMap (type: float)
        lemon::ListGraph::EdgeMap<float> cost(g);
        int index = 0;
        for(int i = 0; i < dist.cols ; i++){
            for(int j = 0; j < dist.rows ; j++){
                int gap = std::abs(gapsizes.at(j));
    
                if( dist.at<float>(j,i) > (gap+1)*maxspeed){
                    cost[g.edgeFromId(index)] = minFloat;
                }
                else{
                    cost[g.edgeFromId(index)] = (float) max + 1.f - dist.at<float>(j,i); //take negative weights to determine the minimum weight matching
                }
                index++;
            }
        }

        /// create object for maximum weighted graph matching
        lemon::MaxWeightedMatching<lemon::ListGraph, lemon::ListGraph::EdgeMap<float> > mwm(g, cost);
        /// run the algorithm
        mwm.run(); 

        std::vector<int> match;
        /// return the matchings

        int id1 = 0, id2 = 0;
        for(int i = 0; i < n_edges; i++){
            lemon::ListGraph::Edge e = g.edgeFromId(i);

            if(mwm.matching(e)){
                if(cost[e] > 0){
                    Tracking::combineTwoCellTracks(tracks1, id_end.at(id1), id_start.at(id2));
                }
            }

            id2++;

            /// increase id1-counter in every xth step (x=#tracks with starting time within the gapsize)
            if(i > 0 && ((i+1) % id_start.size()) == 0){
                id1++;
            }
            else if(i == 0 && id_start.size( ) == 1){
                id1++;
            }

            /// reset id2-counter in every xth step (x=#tracks with starting time within the gapsize)
            if(id2 % id_start.size() == 0){
                id2 = 0;
            }
        }
    }

    void generateCostFunctionPairwiseEuklideanDistanceAndPhagocytosis(cv::Mat &cost, std::vector<RegionP> &regions1, std::vector<RegionP> &regions2){
        
        int K1 = regions1.size();
        int K2 = regions2.size();

        cv::Mat c(K1, K2, CV_32FC1, cv::Scalar(0));
        cv::Mat frs(K1, K2, CV_32FC1, cv::Scalar(0));

        std::vector<RegionP>::iterator regit1, regit2;

        int i = 0;
        for(regit1 = regions1.begin(); regit1 != regions1.end(); regit1++){
            
            int j = 0;
            double fr1 = regit1->getFungalRatioGreen();
            
            for(regit2 = regions2.begin(); regit2 != regions2.end(); regit2++){
                double fr2 = regit2->getFungalRatioGreen();
                double fr;

                fr = (std::abs(fr1-fr2));
    
                frs.at<float>(i,j) = (float) fr;
                c.at<float>(i,j) = (float) computeDistance(*regit1, *regit2);
                j++;
            }
            i++;
        }

        cv::multiply(frs, c, frs);
        cv::add(frs, c, frs);

        frs.copyTo(cost);

    }

    void combineCellTracks(std::vector<CellTrackP> &tracks1, const int t_max, const int maxgapsize){

        std::vector<CellTrackP>::iterator trackit;
        double maxspeed = Tracking::getMaxSpeed(tracks1);
    
        for(int t = t_max; t >= 0 ; t--){

            std::vector< int > id_start, id_start2;
            std::vector< int > id_end, id_end2;
            std::vector< int > gapsizes;
            std::vector<RegionP> reg_start, reg_end, reg_start2, reg_end2;

            for(trackit =  tracks1.begin(); trackit != tracks1.end(); trackit++){

                int t_end = trackit->getEndTime();
                int t_start = trackit->getStartTime();

                if(t_end == t){
                    id_end.push_back(trackit->getID());
                    RegionP tmp;
    
                    trackit->getLastRegion(tmp);
                    reg_end.push_back(tmp);
                }

                if(t_start >= t+1 && t_start <= t+maxgapsize+1){
                    id_start.push_back(trackit->getID());
                    RegionP tmp;
    
                    trackit->getFirstRegion(tmp);
                    reg_start.push_back(tmp);

                    int d = t_start-t-1;
                    gapsizes.push_back(d);
                }
            }

            if(id_start.size() > 0 && id_end.size() > 0){
                Tracking::combineTracksWithGraph(tracks1, id_start, id_end, reg_start, reg_end,  gapsizes, maxspeed);
            }
        }
    }

    void combineCellTracks2(std::vector<CellTrackP> &tracks1, const int t_max, const int maxgapsize){

        std::vector<CellTrackP>::iterator trackit;
        double maxspeed = Tracking::getMaxSpeed(tracks1);
    
        for(int t = t_max; t >= 0 ; t--){

            std::vector< int > id_start, id_start2;
            std::vector< int > id_end, id_end2;
            std::vector< int > gapsizes;
            std::vector<RegionP> reg_start, reg_end, reg_start2, reg_end2;

            for(trackit =  tracks1.begin(); trackit != tracks1.end(); trackit++){

                int t_end = trackit->getEndTime();
                int t_start = trackit->getStartTime();

                if(t_end == t){
                    id_end.push_back(trackit->getID());
                    RegionP tmp;
    
                    trackit->getLastRegion(tmp);
                    reg_end.push_back(tmp);
                }

                if(t_start >= t+1 && t_start <= t+maxgapsize+1){
                    id_start.push_back(trackit->getID());
                    RegionP tmp;

                    trackit->getFirstRegion(tmp);
                    reg_start.push_back(tmp);

                    int d = t_start-t-1;
                    gapsizes.push_back(d);
                }

            }

            if(id_start.size() > 0 && id_end.size() > 0){
                combineTracksWithGraph2(tracks1, id_start, id_end, reg_start, reg_end,  gapsizes, maxspeed);
            }

        }

    }

    void combineTracksWithGraph(std::vector<CellTrack> &tracks1, std::vector<int> &id_start, std::vector<int> &id_end, std::vector<Region> &reg_start, std::vector<Region> &reg_end, std::vector<int> &gapsizes, const double &maxspeed){
    
        cv::Mat dist(id_start.size(), id_end.size(), CV_32F, cv::Scalar(0));
        Tracking::computeDistances(dist, id_start, id_end, reg_start, reg_end);

        // solve with graph theory and declare graph object
        lemon::ListGraph g;

        // add nodes to the graph
        int n_nodes = id_end.size() + id_start.size();
        int n_edges = id_end.size() * id_start.size();

        for(size_t i = 0; i <  n_nodes; i++){
            g.addNode();
        }

        for(size_t i = 0; i < id_end.size(); i++){
            for(size_t j = id_end.size(); j < n_nodes; j++){
                g.addEdge(g.nodeFromId(i), g.nodeFromId(j));
            }
        }

        double min, max;
        cv::minMaxLoc(dist, &min, &max);

        // add weights to an EdgeMap (type: float)
        lemon::ListGraph::EdgeMap<float> cost(g);
        int index = 0;
        for(size_t i = 0; i < dist.cols ; i++){
            for(size_t j = 0; j < dist.rows ; j++){
                int gap = std::abs(gapsizes.at(j));
    
                if( dist.at<float>(j,i) > (gap+1)*maxspeed){
                    cost[g.edgeFromId(index)] = minFloat;
                }
                else{
                    // take negative weights to determine the minimum weight matching
                    cost[g.edgeFromId(index)] = (float) max + 1.f - dist.at<float>(j,i); 
                }
    
                index++;
            }
        }

        // create object for maximum weighted graph matching
        lemon::MaxWeightedMatching<lemon::ListGraph, lemon::ListGraph::EdgeMap<float> > mwm(g, cost);
        // run the algorithm
        mwm.run(); 
        
        std::vector<int> match;
        // return the matchings

        int id1 = 0, id2 = 0;
        for(size_t i = 0; i < n_edges; i++){
            lemon::ListGraph::Edge e = g.edgeFromId(i);

            if(mwm.matching(e)){
                if(cost[e] > 0){
                    Tracking::combineTwoCellTracks(tracks1, id_end.at(id1), id_start.at(id2));
                }
            }

            id2++;

            // increase id1-counter in every xth step (x=#tracks with starting time within the gapsize)
            if(i > 0 && ((i+1) % id_start.size()) == 0){
                id1++;
            } else if(i == 0 && id_start.size( ) == 1){
                id1++;
            }

            // reset id2-counter in every xth step (x=#tracks with starting time within the gapsize)
            if(id2 % id_start.size() == 0){
                id2 = 0;
            }
        }
    }

    /**
     * This function merges two cell tracks to one.
     * The resulting track has the ID of the first cell track and the second track is deleted.
     *
     * @param tracks1 set of cell tracks
     * @param id1 ID of first cell track
     * @param id2 ID of second cell track
     */
    void combineTwoCellTracks(std::vector<CellTrack> &tracks1, const int &id1, const int &id2){

        std::vector<CellTrack>::iterator it1, it2, it3;

        for(it1 = tracks1.begin(); it1 != tracks1.end(); it1++){

            if(it1->getID() == id1){

                for(it2 = tracks1.begin(); it2 != tracks1.end(); it2++){
                    if(it2->getID() == id2){

                        // combine the tracks
                        it1->add(*it2);
                        tracks1.erase(it2);

                        // change cell id for interactions
                        for(it3 = tracks1.begin(); it3 != tracks1.end(); it3++){
                            it3->changeInteractionID(id2, id1);
                        }

                        break;
                    }
                }

                break;
            }
        }

    }

    void combineTwoCellTracks(std::vector<CellTrackP> &tracks1, int id1, int id2){

        std::vector<CellTrackP>::iterator it1, it2, it3;

        for(it1 = tracks1.begin(); it1 != tracks1.end(); it1++){

            if(it1->getID() == id1){

                for(it2 = tracks1.begin(); it2 != tracks1.end(); it2++){
                    if(it2->getID() == id2){

                        /// combine the tracks
                        it1->add(*it2);
                        tracks1.erase(it2);

                        /// change cell id for interactions
                        for(it3 = tracks1.begin(); it3 != tracks1.end(); it3++){
                            it3->changeInteractionID(id2, id1);
                        }
                        break;
                    }

                }

                break;
            }

        }

    }

    void computeDistances(cv::Mat & dist, const std::vector<int> &id_start, const std::vector<int> &id_end, std::vector<Region> &reg_start, std::vector<Region> &reg_end){

        for(size_t i = 0; i < id_start.size(); i++){
            for(size_t j = 0; j < id_end.size(); j++){
                dist.at<float>(i,j) = (float) computeDistance(reg_start.at(i), reg_end.at(j));

            }
        }

    }

    void computeDistances(cv::Mat & dist, const std::vector<int> &id_start, const std::vector<int> &id_end, std::vector<RegionP> &reg_start, std::vector<RegionP> &reg_end){

        for(size_t i = 0; i < id_start.size(); i++){
            for(size_t j = 0; j < id_end.size(); j++){
                dist.at<float>(i,j) = (float) computeDistance(reg_start.at(i), reg_end.at(j));
            }
        }

    }

    void computeDistancesToImageBoundaries(cv::Mat &dist, const std::vector<int> &id_start, const std::vector<int> &id_end, std::vector<Region> &reg_start, std::vector<Region> &reg_end, int rows, int cols){

        for(size_t i = 0; i < id_start.size(); i++){
            for(size_t j = 0; j < id_end.size(); j++){
                dist.at<float>(i,j) = (float) Tracking::computeDistanceToImageBoundaries(reg_start.at(i), reg_end.at(j), rows, cols);
            }
        }

    }

    double computeDistanceToImageBoundaries(Region &r, Region &p, int rows, int cols){
        
        cv::Point center1 = r.getCentroidPoint();

        std::vector<int> elements;

        elements.push_back(center1.x);
        elements.push_back(center1.y);
        elements.push_back(cols - center1.x);
        elements.push_back(rows - center1.y);

        cv::Point center2 = p.getCentroidPoint();
        elements.push_back(center2.x);
        elements.push_back(center2.y);
        elements.push_back(cols - center2.x);
        elements.push_back(rows - center2.y);

        std::vector<int>::iterator it = min_element(elements.begin(), elements.end());
    
        return *it;
    }

    void combineCellTracks(std::vector<CellTrack> &tracks1, const int &t_max, const int &maxgapsize){
        std::vector<CellTrack>::iterator trackit;

        double maxspeed = Tracking::getMaxSpeed(tracks1);

        std::cout << "\t[5.2.1] Max speed during iteration: " << maxspeed << std::endl;
        std::cout << "\t[5.2.2] Max gap size: " << maxgapsize << std::endl;
        std::cout << "\t[5.2.3] MinFloat: " << minFloat << std::endl;
        
        for(int t = t_max; t >= 0 ; t--){

            std::vector<int> id_start, id_start2;
            std::vector<int> id_end, id_end2;
            std::vector<int> gapsizes;
            std::vector<Region> reg_start, reg_end, reg_start2, reg_end2;

            for(trackit =  tracks1.begin(); trackit != tracks1.end(); trackit++){

                int t_end = trackit->getEndTime();
                int t_start = trackit->getStartTime();

                if(t_end == t){
                    id_end.push_back(trackit->getID());
                    Region tmp;
                    trackit->getLastRegion(tmp);
                    reg_end.push_back(tmp);
                }

                if(t_start >= t+1 && t_start <= t+maxgapsize+1){
                    id_start.push_back(trackit->getID());
                    Region tmp;
                    trackit->getFirstRegion(tmp);
                    reg_start.push_back(tmp);

                    int d = t_start-t-1;
                    gapsizes.push_back(d);
                }
            }

            if(id_start.size() > 0 && id_end.size() > 0){
                combineTracksWithGraph(tracks1, id_start, id_end, reg_start, reg_end,  gapsizes, maxspeed);
            }
        }
    }

    void combineCellTracks(std::vector<CellTrack> &tracks1, const int t_max, const int maxgapsize, double maxspeed, int rows, int cols, const bool &DEBUG){

        std::vector<CellTrack>::iterator trackit;

        if(maxspeed == 0){
            maxspeed = Tracking::getMaxSpeed(tracks1);
        }

        for(int t = t_max; t >= 0 ; t--){
    
            std::vector< int > id_start, id_start2;
            std::vector< int > id_end, id_end2;
            std::vector< int > gapsizes;
            std::vector<Region> reg_start, reg_end, reg_start2, reg_end2;

            for(trackit =  tracks1.begin(); trackit != tracks1.end(); trackit++){

                int t_end = trackit->getEndTime();
                int t_start = trackit->getStartTime();

                if(t_end == t){
                    id_end.push_back(trackit->getID());
                    Region tmp;
                    trackit->getLastRegion(tmp);
                    reg_end.push_back(tmp);
                }

                if(t_start >= t+1 && t_start <= t+maxgapsize+1){
                    id_start.push_back(trackit->getID());
                    Region tmp;

                    trackit->getFirstRegion(tmp);
                    reg_start.push_back(tmp);

                    int d = t_start-t-1;
                    gapsizes.push_back(d);
                }

            }

            if(id_start.size() > 0 && id_end.size() > 0){
                combineTracksWithGraph(tracks1, id_start, id_end, reg_start, reg_end,  gapsizes, maxspeed, rows, cols, DEBUG);
            }

        }

    }

    /**
     * In the post-processing step all regions assigned to cluster-class that could not be split in previous steps are split
     * dependent on the number of overlapping singlecells-class regions.
     *
     * @param regions two-dimensional set of regions
     * @param tracks_singlecells set of cell tracks
     * @param tracks_cluster set of cell tracks of class cluster
     */
    void postProcessingClusterSplitting(std::vector<std::vector<Region>> &regions, std::vector<CellTrack> &tracks_singlecells, std::vector<CellTrack> &tracks_cluster){
        std::vector< CellTrack>::iterator track_it, track2_it;
        
        int t_max = regions.size();

        std::vector<CellTrack> track_tmp2;

        /// iterate over all cluster - tracks
        for(track2_it =  tracks_cluster.begin(); track2_it != tracks_cluster.end(); track2_it++){

            Tracking::computeOverlapRegionsForCellTracks(tracks_singlecells, tracks_cluster, t_max);

            int t = track2_it->getStartTime();

            std::vector<Region> overlapb, overlapf;

            int b1 = track2_it->overlapregionsb1.size();
    //		int b2 = track2_it->overlapregionsb2.size();
            int f1 = track2_it->overlapregionsf1.size();
    //		int f2 = track2_it->overlapregionsf2.size();


            if(b1 > 0 || f1 > 0){
                std::vector<Region> overlap;

                if(track2_it->overlapregionsb1.size() > track2_it->overlapregionsf1.size()){
                    overlap = track2_it->overlapregionsb1;
                }
                else{
                    overlap = track2_it->overlapregionsf1;
                }

                Region cluster;
                track2_it->getRegionAt(t, cluster);

                RegionSplitting::separateClusterGMM(cluster, overlap);

                // add regions to cell tracks
                for(size_t i = 0; i < overlap.size(); i++){
                    
                    Tracking::addRegionToTrack(overlap.at(i), tracks_singlecells, overlap.at(i).getId(), t);
                    
                    regions.at(t).push_back(overlap.at(i));
                    
                    if(overlap.size() > 1){
                        Tracking::addInteractionToTrack(tracks_singlecells, overlap.at(i).getId(), t, overlap);
                    }
                }

                // remove cluster region from regions vector
                Tracking::deleteRegion(regions, t, track2_it->getID());
            }
            else{
                track_tmp2.push_back(*track2_it);
            }
        }

        tracks_cluster = track_tmp2;
    }

    /**
     * This function comines cell tracks using matching in bipartite graphs.
     *
     * Iteration over all timepoints specified by t_min and t_max.
     * Then for a time point t all ending tracks are extracted and all starting tracks within t+1 and t+maxgapsize.
     *
     * @param tracks1 set of cell tracks
     * @param t_min minimum time point
     * @param t_max maximum time point
     * @param maxgapsize maximum gap size
     */
    void combineCellTracks(std::vector<CellTrack> &tracks1, const int &t_min, const int &t_max, const int &maxgapsize){

        std::vector<CellTrack>::iterator trackit;
        double maxspeed = Tracking::getMaxSpeed(tracks1);

        for(size_t t = t_max; t >= t_min ; t--){

            std::vector<int> id_start, id_start2;
            std::vector<int> id_end, id_end2;
            std::vector<int> gapsizes;
            std::vector<Region> reg_start, reg_end, reg_start2, reg_end2;

            for(trackit =  tracks1.begin(); trackit != tracks1.end(); trackit++){

                int t_end = trackit->getEndTime();
                int t_start = trackit->getStartTime();

                if(t_end == t){
                    id_end.push_back(trackit->getID());
                    Region tmp;

                    trackit->getLastRegion(tmp);
                    reg_end.push_back(tmp);
                }
   
                if(t_start >= t+1 && t_start <= t+maxgapsize+1){
                    id_start.push_back(trackit->getID());
                    Region tmp;
   
                    trackit->getFirstRegion(tmp);
                    reg_start.push_back(tmp);

                    int d = t_start-t-1;
                    gapsizes.push_back(d);
                }

            }

            if(id_start.size() > 0 && id_end.size() > 0){
                Tracking::combineTracksWithGraph(tracks1, id_start, id_end, reg_start, reg_end,  gapsizes, maxspeed);
            }

        }
    }

    /**
     * this function removes all cell tracks that have a length less or equal than maxgapsize
     */
    void postProcessing(std::vector<CellTrack> &tracks, const int &maxgapsize){

        std::vector<CellTrack>::iterator trackit;
        
        std::vector<CellTrack> tmp;

        for(trackit =  tracks.begin(); trackit != tracks.end(); trackit++){
            int l = trackit->getLength();
        
            int error = trackit->getNumberOfMissingMeasurements();
            int n = l - error;
        
            if(n > maxgapsize){
                tmp.push_back(*trackit);
            }

        }

        tracks = tmp;
    }

    /**
     * this function removes all cell tracks that have a length less or equal than maxgapsize
     */
    void postProcessing(std::vector<CellTrackP> &tracks, int maxgapsize){

        std::vector<CellTrackP>::iterator trackit;
        std::vector<CellTrackP> tmp;

        for(trackit =  tracks.begin(); trackit != tracks.end(); trackit++){
            
            int l = trackit->getLength();
            int error = trackit->getNumberOfMissingMeasurements();
            int n = l - error;
            
            if(n > maxgapsize){
                tmp.push_back(*trackit);
            }

        }

        tracks = tmp;
    }

    /**
     * This function renumbers the IDs of a set of cell tracks and puts them in increasing order without missing values.
     *
     * @param tracks set of cell tracks
     */
    void correctIds(std::vector<CellTrack> &tracks){
        std::vector<CellTrack>::iterator it;
        
        int id = 1;

        for(it = tracks.begin(); it != tracks.end(); it++, id++){
            if(it->getID() != id){
                it->setId(id);
            }
        }

    }

    /**
     * @param tracks set of cell tracks
     * @return number of interactions
     */
    int getNumberOfInteractions(std::vector<CellTrack> &tracks){
        int n = 0;

        std::vector<CellTrack>::iterator trackit = tracks.begin();

        for(; trackit != tracks.end(); trackit++){
            n += trackit->getNumberOfInteractions();
        }

        return n;
    }

} // Tracking