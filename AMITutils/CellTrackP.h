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

#pragma once

#include <string>
#include <vector>
#include "CellTrack.h"
#include "RegionP.h"


class CellTrackP: public CellTrack {
    
    public:
        
        CellTrackP(int t0, int id);
        virtual ~CellTrackP();

        std::vector<RegionP> overlapregionsf1;  /**< vectors containing overlapping regions of class 1 in the next timestep */
        std::vector<RegionP> overlapregionsf2;  /**< vectors containing overlapping regions of class 2 in the next timestep */

        std::vector<RegionP> overlapregionsb1; /**< vectors containing overlapping regions of class 1 in the previous timestep */
        std::vector<RegionP> overlapregionsb2; /**< vectors containing overlapping regions of class 2 in the previous timestep */

        std::vector<RegionP> overlapregionsf1_fc;  /**< vectors containing overlapping regions of class 1 in the next timestep */
        std::vector<RegionP> overlapregionsf2_fc;  /**< vectors containing overlapping regions of class 2 in the next timestep */

        std::vector<RegionP> overlapregionsb1_fc; /**< vectors containing overlapping regions of class 1 in the previous timestep */
        std::vector<RegionP> overlapregionsb2_fc; /**< vectors containing overlapping regions of class 2 in the previous timestep */

        void add(RegionP &r);
        void add(std::vector< RegionP > & regions);
        void add(CellTrackP &track1);
        void add(RegionP &region, int t);
        void addFungiId(int t, int id);
        void addPhagocytosis(int tstart, int id);
        void addPhagocytosis(int tstart, int id, int pos);

        void change(RegionP &r, int t);
        void changeFID(int id, int id_new);
        void clearOverlaps();

        void deleteFirstRegion();
        void deleteLastRegion();
        void deleteRegion(int t);

        bool existsAt(int t);

        int getAverageArea();
        const int getEndTime();
        void getFIDs(int t, std::vector<int> &fids);
        int getFirstPhagocytosisTimepoint();
        void getFirstRegion(RegionP &region);
        int getIState(int t);
        void getIStates(std::vector<int> &states);
        void getLastRegion(RegionP &region);
        int getLength();
        double getMaxSpeed();
        int getNumberOfFungi(int t);
        int getNumberOfIState(int i);
        int getNumberOfMissingMeasurements();
        int getNumberOfPhagocytoses(int t);
        int getNumberOfPhagocytozedFungi(int t);
        void getPhagocytosisIds(int t, std::vector<int> &ids);
        RegionP * getRegionAt(int t);
        void getRegionAt(int t, RegionP & region);
        void getTimepointsOfPhagocytoses(std::vector<int> & timepoints);
        void getTrack(std::vector<RegionP> &track);

        bool hasFungiId(int id);
        bool hasFungiId(int id, int t);
        bool hasPhagocytozedFungus(int id);
        bool hasValidMeasurement(int t);

        bool isInteracting(int t);
        bool isTouchingMoreFungiThan(int t, int fid);

        void printToFile(std::string path);

        void removeFungiIds(int t);
        void removeFungiId(int t, int id);
        void removeLastPhID(int t);
        void removeNaNRegions();
        void removeRegionAt(int t);
        void removeRegions(int t);

        void setInteraction(int t, int id);
        void setIState(int t, int state);

        /**
         * A cell track is defined to be < than another cell track, if its ending time is less than the other cell track.
         */
        friend bool operator < ( CellTrackP track1,  CellTrackP track2){
            int tend1 = track1.getEndTime();
            int tend2  = track2.getEndTime();
            return tend1 < tend2;
        }

    private:
        std::vector<RegionP> track;

};