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

 #include <vector>
 #include "Region.h"


 class CellTrack
 {
    private:
        std::vector<Region> track;          /* vector containing measurements of the track */

    protected:
        int ID;                             /* CellTrack ID */
	    int t0;                             /* starting time */

    public:

        std::vector<Region> overlapregionsf1;  /**< vectors containing overlapping regions of class singlecells in the next timestep */
        std::vector<Region> overlapregionsf2;  /**< vectors containing overlapping regions of class cluster in the next timestep */

        std::vector<Region> overlapregionsb1; /**< vectors containing overlapping regions of class singlecells in the previous timestep */
        std::vector<Region> overlapregionsb2; /**< vectors containing overlapping regions of class cluster in the previous timestep */

        std::vector<Region> overlapregionsf1_fc;  /**< vectors containing overlapping regions of class singlecells in the next timestep */
        std::vector<Region> overlapregionsf2_fc;  /**< vectors containing overlapping regions of class cluster in the next timestep */

        std::vector<Region> overlapregionsb1_fc; /**< vectors containing overlapping regions of class singlecells in the previous timestep */
        std::vector<Region> overlapregionsb2_fc; /**< vectors containing overlapping regions of class cluster in the previous timestep */

        // constructor
        CellTrack(int t0, int id);
        // destructor
	    virtual ~CellTrack();

        // class function members 
        virtual void add(Region &r);        
        virtual void add(Region &region, const int &t);
        virtual void add(std::vector<Region> &regions);
        virtual void add(CellTrack &track1);
        virtual void change(const Region &region, const int &t);

        bool hasValidMeasurement(const int &t);
        bool existsAt(const int &t); // true if track exists at timepoint t, else false
        
        int getID();
        int getLength();
        int getStartTime();
        const int getEndTime();
        void getRegionAt(const int &t, Region &region);
        Region* getRegionAt(const int &t);
        int getAverageArea();
        double getMaxSpeed();
        void getFirstRegion(Region &region);
        void getLastRegion(Region &region);
        void getTrack(std::vector<Region> &track);
        int getNumberOfMissingMeasurements();
        int getNumberOfInteractions();
        bool getInteraction(const int &t);
        void getPrincipalAxesMinEnclosingEllipse(int t, cv::Size2f &size);

        void setInteraction(int t, int id);
        void setId(int id);
        void setStartTime(int t);

        void clearOverlaps();
        void changeInteractionID(int id, int id_new);

        void removeRegions(int t);

        void printToFile(const std::string path);
        
 };