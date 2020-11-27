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

#include "Outputs.h"
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include "IOputs.h"


namespace outputs 
{

    /**
     * This function prints fungi tracks to one file.
     *
     * @param path path to ouput file
     * @param tracks set of fungi tracks
     */
    void printDataSetFungi(std::string path, std::vector<CellTrack> &tracks){

        std::vector<CellTrack>::iterator trackit;

        ofstream outfile;
        outfile.open(path.c_str());

        if(outfile.is_open()){
            outfile << "t\tx\ty\tID\tA\tn_regions" << std::endl;
        }

        for(trackit = tracks.begin(); trackit != tracks.end(); trackit++){
            /// get initial time
            int t = trackit->getStartTime(); 
            std::vector<Region> regions;
            std::vector<Region>::iterator regit;
            
            trackit->getTrack(regions);
            int ID = trackit->getID();

            for(regit = regions.begin(); regit != regions.end(); regit++,t++){
                cv::Scalar c = regit->getCentroidScalar();
    
                if(outfile.is_open()){
                    /// increase time for the next step
                    outfile << t << "\t"; 
                    outfile << c[1] << "\t";
                    outfile << c[0] << "\t";
                    outfile << ID << "\t";
                    outfile << regit->getArea();
                    outfile << "\t";

                    if(!trackit->hasValidMeasurement(t)){
                        outfile << "nan" ;
                    }
                    else{
                        outfile << regit->n_regions;
                    }
                    outfile << std::endl;
                }
            }
        }

        outfile.close();

    }

    /**
     * This function prints the cell tracks including information about phagocytosis events to one file.
     * The header of the file contains the following parameters:
     * t: time \n
     * x: x-coordinate \n
     * y: y-coordinate \n
     * ID \n
     * A: area \n
     * state: [0-5] state of the cell \n
     * nPh: number of phagocytoses until the current time point \n
     * interaction: [0/1] boolean value for cell interacting with another cell \n
     * interaction IDs: list of IDs of interacting cells \n
     * FIDs: list of fungi IDs that are touching/phagocytozed \n
     * phIDs: list of fungi IDs that are phagocytozed
     *
     * @param path path to outputfile
     * @param tracks set of cell tracks
     */
    void printDataSetPhagocytosis(std::string path, std::vector<CellTrackP> &tracks){

        std::vector<CellTrackP>::iterator trackit;

        ofstream outfile;
        outfile.open(path.c_str());

        if(outfile.is_open()){
            outfile << "t\tx\ty\tID\tA\tstate\tnPh\tinteraction\tinteractionIDs\tFIDs\tdFIDs\tphIDs" << std::endl;
        }

        for(trackit = tracks.begin(); trackit != tracks.end(); trackit++){
            /// get initial time
            int t = trackit->getStartTime(); 
            std::vector<RegionP> regions;
            std::vector<RegionP>::iterator regit;
            
            trackit->getTrack(regions);
            int ID = trackit->getID();

            for(regit = regions.begin(); regit != regions.end(); regit++){
                
                cv::Scalar c = regit->getCentroidScalar();
                int state = regit->getIState();

                if(outfile.is_open()){
                    /// increase time for the next step
                    outfile << t++ << "\t"; 
                    outfile << c[1] << "\t";
                    outfile << c[0] << "\t";
                    outfile << ID << "\t";
                    outfile << trackit->getAverageArea() << "\t";
                    outfile << state << "\t";
                    outfile << regit->getNumberOfPhagocytoses();

                    if(regit->isInteracting()){
                        outfile << "\t" << 1;
                    }
                    else{
                        outfile << "\t" << 0;
                    }

                    /// print interacting PMN IDs
                    std::vector<int> ids;
                    regit->getInteractionIds(ids);
                    if(ids.size() == 0){
                        outfile << "\t" << 0;
                    }
                    else{
                        outfile << "\t";
                        for(size_t i = 0; i < ids.size()-1; i++){
                            outfile << ids.at(i) << ",";
                        }
                        outfile << ids.at(ids.size()-1);
                    }

                    /// print interacting fungi IDs
                    ids.resize(0);
                    regit->getFungiIds(ids);
                    if(ids.size() == 0){
                        outfile << "\t" << 0;
                    }
                    else{
                        outfile << "\t";
                        for(size_t i = 0; i < ids.size()-1; i++){
                            outfile << ids.at(i) << ",";
                        }
                        outfile << ids.at(ids.size()-1);
                    }

                    /// print dead fungi IDs
                    ids.resize(0);
                    regit->getDeadFungiIds(ids);
                    if(ids.size() == 0){
                        outfile << "\t" << 0;
                    }
                    else{
                        outfile << "\t";
                        for(size_t i = 0; i < ids.size()-1; i++){
                            outfile << ids.at(i) << ",";
                        }
                        outfile << ids.at(ids.size()-1);
                    }

                    /// print phagcoytozed fungi IDs
                    ids.resize(0);
                    regit->getPhagocytosisIds(ids);
                    if(ids.size() == 0){
                        outfile << "\t" << 0;
                    }
                    else{
                        outfile << "\t";
                        for(size_t i = 0; i < ids.size()-1; i++){
                            outfile << ids.at(i) << ",";
                        }
                        outfile << ids.at(ids.size()-1);
                    }
                    outfile << std::endl;
                }
            }
        }

        outfile.close();

    }

    void showSegmentedRegions2D(std::string &outputdir, std::string &dir_orig_images , int delay, std::vector<std::vector<Region*>> &regions, bool idoutput) {
    
        std::vector<cv::Vec2i>::iterator vecit;
	    std::vector<Region*>::iterator regit;

        (void) io::create_directory(outputdir);
        
        std::vector<std::string> files_orig;
        io::read_filenames(dir_orig_images, files_orig);
    
        int max = std::min(files_orig.size(), regions.size()); 

        for(size_t i = 0; i < max; i++){
            
            std::string file = files_orig.at(i+delay);
            
            //read as color image
            cv::Mat image = cv::imread(file, 1); 

            /// draw regions
            for(regit = regions.at(i).begin(); regit != regions.at(i).end(); regit++){

                /// contour
                for(vecit = (*regit)->contour_pixels.begin(); vecit != (*regit)->contour_pixels.end(); vecit++){
                    int  i = vecit->val[0];
                    int  j = vecit->val[1];
                    
                    cell_class::immune klass = (*regit)->getClass();

                    if(klass == cell_class::immune::CLUSTER ){
                        image.at<cv::Vec3b>(i,j)[0] = 0; //blue
                        image.at<cv::Vec3b>(i,j)[1] = 0; //green
                        image.at<cv::Vec3b>(i,j)[2] = 255; //red
                    }
                    else if(klass == cell_class::immune::NOISE ){
                        image.at<cv::Vec3b>(i,j)[0] = 255; //blue
                        image.at<cv::Vec3b>(i,j)[1] = 255; //green
                        image.at<cv::Vec3b>(i,j)[2] = 255; //red
                    }
                    else{
                        image.at<cv::Vec3b>(i,j)[0] = 255; //blue
                        image.at<cv::Vec3b>(i,j)[1] = 0; //green
                        image.at<cv::Vec3b>(i,j)[2] = 0; //red
                    }
                }

                /// centroid
                cv::Point center = (*regit)->getCentroidPoint();
                cv::circle(image, center, 1, cv::Scalar(255,255,255), -1, 8, 0 );

                /// ID
                if(idoutput){
                    int id = (*regit)->getId();
                    std::stringstream ss;
                    ss << id;
                    cv::putText(image, ss.str(), center, cv::FONT_HERSHEY_SIMPLEX, 1, cv::Scalar(255,255,255));
                }
            }

            // save the images
            io::write_image(outputdir, image, (i+delay), false, std::nullopt, "jpg" );

        }
    }

    void showClassifiedSegmentedRegions2D(std::string &outputdir, const std::string &dir_orig_images , const int &delay, std::vector<std::vector<Region>> &regions, const bool &classoutput, const bool &idoutput){

        std::vector<cv::Vec2i>::iterator vecit;
        std::vector<Region>::iterator regit;

        (void) io::create_directory(outputdir);

		// extract filename of specified directory
        std::vector<std::string> files_orig;
		
		path p(dir_orig_images);
		if(is_directory(p)) {
			for(const auto& entry : boost::make_iterator_range(directory_iterator(p), {})) {	
				files_orig.push_back(entry.path().string());
            }
		}   
        // take care, that images are sorted!
		std::sort(files_orig.begin(), files_orig.end());

        size_t max = std::min(files_orig.size(), regions.size());

        for(size_t i = 0; i < max; i++){
            std::string file = files_orig.at(i+delay);
            // read as color image
            cv::Mat image = cv::imread(file, 1); 

            // draw regions
            for(regit = regions.at(i).begin(); regit != regions.at(i).end(); regit++){

                for(vecit = regit->contour_pixels.begin(); vecit != regit->contour_pixels.end(); vecit++){
                    int i = vecit->val[0];
                    int j = vecit->val[1];
                    
                    cell_class::immune klass = regit->getClass();

                    if(klass == cell_class::immune::CLUSTER ){
                        image.at<cv::Vec3b>(i,j)[0] = 0;    //blue
                        image.at<cv::Vec3b>(i,j)[1] = 0;    //green
                        image.at<cv::Vec3b>(i,j)[2] = 255;  //red
                    }
                    else if(klass == cell_class::immune::NOISE ){
                        image.at<cv::Vec3b>(i,j)[0] = 255;  //blue
                        image.at<cv::Vec3b>(i,j)[1] = 255;  //green
                        image.at<cv::Vec3b>(i,j)[2] = 255;  //red
                    }
                    else{
                        image.at<cv::Vec3b>(i,j)[0] = 255;  //blue
                        image.at<cv::Vec3b>(i,j)[1] = 0;    //green
                        image.at<cv::Vec3b>(i,j)[2] = 0;    //red
                    }                    
                }

                // centroid
                cv::Point center = regit->getCentroidPoint();
                cv::circle(image, center, 1, cv::Scalar(255,255,255), -1, 8, 0 ); 

                // ID
                std::stringstream ss;
                
                if(idoutput){
                    int id = regit->getId();
                    ss << id;
                }
                if(classoutput){
                    cell_class::immune klass = regit->getClass();
                    if(idoutput){
                        ss << " (" << cell_class::getCellClass(klass) << " )";
                    } 
                    else {
                        ss << cell_class::getCellClass(klass);
                    }
                }
                if(classoutput || idoutput){
                    cv::putText(image, ss.str(), center, cv::FONT_HERSHEY_SIMPLEX, 1, cv::Scalar(255,255,255), 2);
                }
            }

            // save the images ( +1 because start suffix at 1 )
            io::write_image(outputdir, image, i+delay+1, false, std::nullopt, "jpg" );

        }
    }


    /************************************************************/
    //          showSegmentedRegions2DAndTrajectory             //
    /************************************************************/

    /**
     * Plots cell contours and trajectories on original images and saves them to an output directory.
     * [4.0] NNA tracking
     *
     * @param outputdir output file
     * @param dir_orig_images path to original images
     * @param delay time delay
     * @param tracks set of cell tracks
     */
    void showSegmentedRegions2DAndTrajectory(const std::string &outputdir, const std::string &dir_orig_images , const int &delay, int maxt, std::vector<CellTrack> &tracks, const cv::Scalar &trackcolor, const cv::Scalar &idcolor, const bool &showInteractions, const cv::Size2i imagesize){
		
        std::vector<cv::Vec2i>::iterator vecit;
		std::vector<CellTrack>::iterator it_track;

		(void) io::create_directory(outputdir);
		std::vector<std::string> files_orig;
		
        path p(dir_orig_images);
		if(is_directory(p)) {
			for(const auto& entry : boost::make_iterator_range(directory_iterator(p), {}))			
				files_orig.push_back(entry.path().string());
		}
        // take care, that images are sorted!
		std::sort(files_orig.begin(), files_orig.end());	

		if(maxt == 0){
			maxt = files_orig.size() - delay;
		}
		else{
			maxt = std::min((int) files_orig.size(), maxt+1)-delay;
		}

		for(size_t t = 0; t < maxt; t++){

            std::string file = files_orig.at(t+delay);
            			
			// read as color image
            cv::Mat image = cv::imread(file, 1);

			// draw regions
			for(it_track = tracks.begin(); it_track != tracks.end(); it_track++){

				if(it_track->hasValidMeasurement(t)){

					// plot trajectories -> centers of mass of all cells in previous frames
                    outputs::plotTrajectory(image, (*it_track), t, trackcolor);

					Region region;
					it_track->getRegionAt(t, region);
					region.computeCentroid();

					// contour
					for(vecit = region.contour_pixels.begin(); vecit != region.contour_pixels.end(); vecit++){
						int  i = vecit->val[0];
						int  j = vecit->val[1];
						image.at<cv::Vec3b>(i,j)[0] = region.color_border[0]; //blue
						image.at<cv::Vec3b>(i,j)[1] = region.color_border[1]; //green
						image.at<cv::Vec3b>(i,j)[2] = region.color_border[2]; //red
					}

					// centroid
					cv::Point center = region.getCentroidPoint();
					cv::circle(image, center, 1, trackcolor, -1, 8, 0 );

					// ID
					int id = it_track->getID();
					std::stringstream ss;
					ss << id;

					std::string id_text = ss.str();
					if(showInteractions && region.isInteracting()){
						id_text += "*";
					}

					cv::putText(image, id_text, center, cv::FONT_HERSHEY_SIMPLEX, 1, idcolor, 2);
				}
			}
			            
			if(imagesize.area() > 0){
                // resize image to imagesize
				cv::resize(image, image, imagesize); 
			}

            // save the images ( +1 because start suffix at 1 )
            io::write_image(outputdir, image, (t+delay)+1, false, std::nullopt, "jpg" );

		}
    }

    /// [5.0] Interaction tracking  &&  [6.0] Combine tracklets globally
    void showSegmentedRegions2DAndTrajectory(const std::string &outputdir, const std::string &dir_orig_images , const int &delay, std::vector<CellTrack> &tracks_single, std::vector<CellTrack> &tracks_cluster){
        std::vector<cv::Vec2i>::iterator it_vec;
        std::vector<CellTrack>::iterator it_track;

        (void) io::create_directory(outputdir);

        std::vector<std::string> files_orig;
        io::read_filenames(dir_orig_images, files_orig);

        for(int t = 0; t < files_orig.size()-1; t++){

            std::string file = files_orig.at(t+delay);

            // read as color image
            cv::Mat image = cv::imread(file, 1);

            // draw regions
            for(it_track = tracks_single.begin(); it_track != tracks_single.end(); it_track++){

                if(it_track->hasValidMeasurement(t)){

                    Region region;
                    it_track->getRegionAt(t, region);
                    region.computeCentroid();

                    // contour
                    for(int ii = 0; ii < region.contour_pixels.size(); ii++){
                        int i = region.contour_pixels.at(ii)[0];
                        int j = region.contour_pixels.at(ii)[1];
                        image.at<cv::Vec3b>(i,j)[0] = region.color_border[0]; //blue
                        image.at<cv::Vec3b>(i,j)[1] = region.color_border[1]; //green
                        image.at<cv::Vec3b>(i,j)[2] = region.color_border[2]; //red
                    }

                    // centroid
                    cv::Point center = region.getCentroidPoint();
                    cv::circle(image, center, 1, cv::Scalar(255,255,255), -1, 8, 0 );

                    // ID
                    int id = it_track->getID();
                    std::stringstream ss;
                    ss << id;
                    cv::putText(image, ss.str(), center, cv::FONT_HERSHEY_SIMPLEX, 1, cv::Scalar(0,0,255), 2);


                    // plot trajectories -> centers of mass of all cells in previous frames
                    outputs::plotTrajectory(image, (*it_track), t);

                }
            }

            // draw regions of cluster-class
            for(it_track = tracks_cluster.begin(); it_track != tracks_cluster.end(); it_track++){

                if(it_track->existsAt(t)){

                    Region region;
                    it_track->getRegionAt(t, region);
                    region.computeCentroid();

                    // contour
                    for(int ii = 0; ii < region.contour_pixels.size(); ii++){
                        int i = region.contour_pixels.at(ii)[0];
                        int j = region.contour_pixels.at(ii)[1];
                        image.at<cv::Vec3b>(i,j)[0] = 255; //blue
                        image.at<cv::Vec3b>(i,j)[1] = 255; //green
                        image.at<cv::Vec3b>(i,j)[2] = 255; //red
                    }

                    // centroid
                    cv::Point center = region.getCentroidPoint();
                    cv::circle(image, center, 1, cv::Scalar(255,255,255), -1, 8, 0 );

                    // ID
                    int id = it_track->getID();
                    std::stringstream ss;
                    ss << id;
                    cv::putText(image, ss.str(), center, cv::FONT_HERSHEY_SIMPLEX, 1, cv::Scalar(0,0,255), 2);

                    // plot trajectories -> centers of mass of all cells in previous frames
                    outputs::plotTrajectory(image, (*it_track), t);

                }

            }

            // save the images ( +1 because start suffix at 1 )
            io::write_image(outputdir, image, (t+delay)+1, false, std::nullopt, "jpg" );

        }

    }

    void showSegmentedRegions2DAndTrajectory(std::string outputdir, std::string dir_orig_images , int delay,  std::vector<CellTrackP> &tracks, std::vector<CellTrackP> &tracks2){
		
        std::vector<cv::Vec2i>::iterator it_vec;
		std::vector<CellTrackP>::iterator it_track;

        (void) io::create_directory(outputdir);
		
        std::vector<std::string> files_orig;
		io::read_filenames(dir_orig_images, files_orig);
        
		for(size_t t = 0; t < files_orig.size()-1; t++){

			std::string file = files_orig.at(t+delay);
			
            /// read as color image
            cv::Mat image = cv::imread(file, 1); 

			/// draw regions
			for(it_track = tracks.begin(); it_track != tracks.end(); it_track++){

				if(it_track->hasValidMeasurement(t)){

					RegionP region;
					it_track->getRegionAt(t, region);
					region.computeCentroid();

					/// contour
					for(size_t ii = 0; ii < region.contour_pixels.size(); ii++){
						int  i = region.contour_pixels.at(ii)[0];
						int  j = region.contour_pixels.at(ii)[1];
						
                        image.at<cv::Vec3b>(i,j)[0] = region.color_border[0]; //blue
						image.at<cv::Vec3b>(i,j)[1] = region.color_border[1]; //green
						image.at<cv::Vec3b>(i,j)[2] = region.color_border[2]; //red
					}

					/// centroid
					cv::Point center = region.getCentroidPoint(); //(region.centroid[1], region.centroid[0]);
					cv::circle(image, center, 1, cv::Scalar(255,255,255), -1, 8, 0 );

					/// ID
					int id = it_track->getID();
					std::stringstream ss;
					ss << id;
					cv::putText(image, ss.str(), center, cv::FONT_HERSHEY_SIMPLEX, 1, cv::Scalar(255,255,255));

					/// plot trajectories -> centers of mass of all cells in previous frames
					outputs::plotTrajectory(image, (*it_track), t);

				}
			}

			/// draw regions of class 2
			for(it_track = tracks2.begin(); it_track != tracks2.end(); it_track++){

				if(it_track->existsAt(t)){

					RegionP region;
					it_track->getRegionAt(t, region);
					region.computeCentroid();

					/// contour
					for(unsigned int ii = 0; ii < region.contour_pixels.size(); ii++){
						int  i = region.contour_pixels.at(ii)[0];
						int  j = region.contour_pixels.at(ii)[1];
						
                        image.at<cv::Vec3b>(i,j)[0] = 255; //blue
						image.at<cv::Vec3b>(i,j)[1] = 255; //green
						image.at<cv::Vec3b>(i,j)[2] = 255; //red
					}

					/// centroid
					cv::Point center = region.getCentroidPoint(); 
					cv::circle(image, center, 1, cv::Scalar(255,255,255), -1, 8, 0 );

					/// ID
					int id = it_track->getID();
					std::stringstream ss;
					ss << id;
					cv::putText(image, ss.str(), center, cv::FONT_HERSHEY_SIMPLEX, 1, cv::Scalar(255,255,255));

					/// plot trajectories -> centers of mass of all cells in previous frames
					outputs::plotTrajectory(image, (*it_track), t);

				}
			}

            // save the images ( +1 because start suffix at 1 )
            io::write_image(outputdir, image, (t+delay)+1, false, std::nullopt, "jpg" );
		
        }
    }

    /// [8.0] Outputs
    void showSegmentedRegions2DAndTrajectory(std::string &outputdir, std::string dir_orig_images, const int &delay, const int &maxt, std::vector<CellTrack> &tracks, const cv::Scalar &trackcolor, const cv::Scalar &idcolor, const bool showInteractions, cv::Size2i imagesize, const int fontsize){
        std::vector<CellTrack*> trackp;

        for(std::vector<CellTrack>::iterator it = tracks.begin(); it != tracks.end(); it++){
            trackp.push_back(&*it);
        }

        outputs::showSegmentedRegions2DAndTrajectory(outputdir, dir_orig_images, delay, maxt, trackp, trackcolor, idcolor, showInteractions, imagesize, fontsize);

    }

    /// [8.0] Outputs
    void showSegmentedRegions2DAndTrajectory(std::string &outputdir, std::string &dir_orig_images, const int delay, int maxt, std::vector<CellTrack*> &tracks, const cv::Scalar &trackcolor, const cv::Scalar &idcolor, const bool showInteractions, const cv::Size2i imagesize, const int fontsize){
        
        std::vector<cv::Vec2i>::iterator vecit;
        std::vector<CellTrack*>::iterator it_track;

        (void) io::create_directory(outputdir);

        std::vector<std::string> files_orig;

        path p(dir_orig_images);
        if(is_directory(p)) {
            for(const auto& entry : boost::make_iterator_range(directory_iterator(p), {}))
                files_orig.push_back(entry.path().string());
        }
        // take care, that images are sorted!
        std::sort(files_orig.begin(), files_orig.end());

        if(maxt == 0){
            maxt = files_orig.size() - delay;
        }
        else{
            maxt = std::min((int) files_orig.size(), maxt+1)-delay;
        }

        for(int t = 0; t < maxt; t++){

            std::string file = files_orig.at(t+delay);

            // read as color image
            cv::Mat image = cv::imread(file, 1);

            //draw regions
            for(it_track = tracks.begin(); it_track != tracks.end(); it_track++){

                if((*it_track)->hasValidMeasurement(t)){

                    // plot trajectories -> centers of mass of all cells in previous frames
                    outputs::plotTrajectory(image, (*it_track), t, trackcolor);

                    Region * region = (*it_track)->getRegionAt(t);
                    region->computeCentroid();

                    // contour
                    for(vecit = region->contour_pixels.begin(); vecit != region->contour_pixels.end(); vecit++){
                        int  i = vecit->val[0];
                        int  j = vecit->val[1];
                        image.at<cv::Vec3b>(i,j)[0] = region->color_border[0]; //blue
                        image.at<cv::Vec3b>(i,j)[1] = region->color_border[1]; //green
                        image.at<cv::Vec3b>(i,j)[2] = region->color_border[2]; //red
                    }

                    // centroid
                    cv::Point center = region->getCentroidPoint(); //(region.centroid[1], region.centroid[0]);
                    cv::circle(image, center, 1, trackcolor, -1, 8, 0 );

                    // ID
                    int id = (*it_track)->getID();
                    std::stringstream ss;
                    ss << id;

                    std::string id_text = ss.str();
                    if(showInteractions && region->isInteracting()){
                        id_text += "*";
                    }

                    cv::putText(image, id_text, center, cv::FONT_HERSHEY_SIMPLEX,  fontsize, idcolor);
                }
            }

            if(imagesize.area() > 0){
                // resize image to imagesize
                cv::resize(image, image, imagesize);
            }

            // save the images ( +1 because start suffix at 1 )
            io::write_image(outputdir, image, (t+delay)+1, false, std::nullopt, "jpg" );

        }
    }

    /************************************************************/
    //                      plotTrajectory                      //
    /************************************************************/

    void plotTrajectory(cv::Mat &image, CellTrack &track, const int &t, const cv::Scalar &trackcolor){

        int t1 = track.getStartTime();
        int t2 = t1+1;
        int tmax = t;

        outputs::plotTrajectory(image, track, t1, t2, tmax, trackcolor);

    }

    void plotTrajectory(cv::Mat &image, CellTrack *track, const int &t, const cv::Scalar &trackcolor){

        int t1 = track->getStartTime();
        int t2 = t1+1;
        int tmax = t;

        outputs::plotTrajectory(image, track, t1, t2, tmax, trackcolor);

    }

    /**
     * plot trajectorie between two timepoints t1 and t2 with t1 < t2
     */
    void plotTrajectory(cv::Mat &image, CellTrack &track, const int &t1, const int &t2, const int &tmax, const cv::Scalar &trackcolor){

        if(t1 <= tmax && t2 <= tmax && track.hasValidMeasurement(t1)){

            if(!track.hasValidMeasurement(t2)){
                outputs::plotTrajectory(image, track, t1, t2+1, tmax, trackcolor);
            }
            else{
                Region r1, r2;
                track.getRegionAt(t1,r1);
                track.getRegionAt(t2,r2);
                cv::Point center1, center2;
                center1 = r1.getCentroidPoint();
                center2 = r2.getCentroidPoint();

                cv::circle(image, center1, 1, trackcolor, -1, 8, 0);
                cv::circle(image, center2, 1, trackcolor, -1, 8, 0);
                cv::line(image, center1, center2, trackcolor, 1);
                outputs::plotTrajectory(image, track, t2, t2+1, tmax, trackcolor);
            }
        }
    }

    /**
     * plot trajectorie between two timepoints t1 and t2 with t1 < t2
     */
    void plotTrajectory(cv::Mat &image, CellTrack *track, const int &t1, const int &t2, const int &tmax, const cv::Scalar &trackcolor){

        if(CellTrackP* ct = dynamic_cast<CellTrackP*>(track)){
            outputs::plotTrajectory(image, ct, t1, t2, tmax, trackcolor);
        }

        if(t1 <= tmax && t2 <= tmax && track->hasValidMeasurement(t1)){

            if(!track->hasValidMeasurement(t2)){
                outputs::plotTrajectory(image, track, t1, t2+1, tmax, trackcolor);
            }
            else{
                Region * r1, * r2;
                r1 = track->getRegionAt(t1);
                r2 = track->getRegionAt(t2);
                cv::Point center1, center2;
                center1 = r1->getCentroidPoint();
                center2 = r2->getCentroidPoint();

                cv::circle(image, center1, 1, trackcolor, -1, 8, 0);
                cv::circle(image, center2, 1, trackcolor, -1, 8, 0);
                cv::line(image, center1, center2, trackcolor, 1);
                outputs::plotTrajectory(image, track, t2, t2+1, tmax, trackcolor);
            }
        }

    }

    /**
     * plot trajectorie between two timepoints t1 and t2 with t1 < t2
     */
    void plotTrajectory(cv::Mat &image, CellTrackP *track, const int t1, const int t2, const int tmax, const cv::Scalar trackcolor){

        if(t1 <= tmax && t2 <= tmax && track->hasValidMeasurement(t1)){

            if(!track->hasValidMeasurement(t2)){
                plotTrajectory(image, track, t1, t2+1, tmax, trackcolor);
            }
            else{
                Region * r1, * r2;
                r1 = track->getRegionAt(t1);
                r2 = track->getRegionAt(t2);
                
                cv::Point center1, center2;
                center1 = r1->getCentroidPoint();
                center2 = r2->getCentroidPoint();

                cv::circle(image, center1, 1, trackcolor, -1, 8, 0);
                cv::circle(image, center2, 1, trackcolor, -1, 8, 0);
                cv::line(image, center1, center2, trackcolor, 1);
                
                plotTrajectory(image, track, t2, t2+1, tmax, trackcolor);
            }
        }
    }

    void plotTrajectory(cv::Mat &image, CellTrack &track, const int &t1, const int &t2, const int &tmax){

        cv::Scalar trackcolor(255,255,255);

        if(t1 <= tmax && t2 <= tmax && track.hasValidMeasurement(t1)){

            if(!track.hasValidMeasurement(t2)){
                outputs::plotTrajectory(image, track, t1, t2+1, tmax, trackcolor);
            }
            else{
                Region r1, r2;
                track.getRegionAt(t1,r1);
                track.getRegionAt(t2,r2);
                cv::Point center1, center2;
                center1 = r1.getCentroidPoint();
                center2 = r2.getCentroidPoint();

                cv::circle(image, center1, 1, trackcolor, -1, 8, 0);
                cv::circle(image, center2, 1, trackcolor, -1, 8, 0);
                cv::line(image, center1, center2, trackcolor, 1);
                outputs::plotTrajectory(image, track, t1+1, t2+1, tmax, trackcolor);
            }
        }
    }

    void plotTrajectory(cv::Mat &image, CellTrack &track, const int &t){

        int t1 = track.getStartTime();
        int t2 = t1+1;
        int tmax = t;

        outputs::plotTrajectory(image, track, t1, t2, tmax);

    }

    void plotAllTrajectories(const std::string &outfile, std::vector<CellTrack> &tracks, const int &cols, const int &rows){

        cv::Mat image(rows, cols, CV_8UC3, cv::Scalar(255,255,255));

        // blue
        cv::Scalar color(255,0,0);

        std::vector<CellTrack>::iterator ct_iterator;
        std::vector<Region>::iterator regit;

        for(ct_iterator = tracks.begin(); ct_iterator != tracks.end(); ct_iterator++){
            std::vector<Region> temp;
            ct_iterator->getTrack(temp);

            cv::Scalar centroid = temp.at(0).getCentroidScalar();

            int x0 = (int) centroid[0];
            int y0 = (int) centroid[1];
            cv::Point p(y0,x0);
            color = cv::Scalar(rand() % 255, rand() % 255, rand() % 255);

            cv::circle(image, p, 1, color, -1, 8, 0);

            for(regit = temp.begin()+1; regit != temp.end(); regit++){
                cv::Scalar c = regit->getCentroidScalar();
                int x = (int) c[0];
                int y = (int) c[1];

                if(!isnan((double)(c[0])) && !isnan((double) c[1])){
                    cv::Point q(y,x);

                    cv::circle(image, q, 1, color, -1, 8, 0);
                    cv::line(image, p, q, color, 1);

                    p = q;
                    x0 = x;
                    y0 = y;
                }
            }
        }

        // save image to specified path and filename
        cv::imwrite(outfile, image);

    }

    void plotAllTrajectories(const std::string &outfile, std::vector<CellTrackP> &tracks, int cols, int rows){

        cv::Mat image(rows, cols, CV_8UC3, cv::Scalar(255,255,255));
        cv::Scalar color(255,0,0); // blue
        
        std::vector<CellTrackP>::iterator ct_iterator;
        std::vector<RegionP>::iterator regit;

        for(ct_iterator = tracks.begin(); ct_iterator != tracks.end(); ct_iterator++){
            
            std::vector<RegionP> temp;
            ct_iterator->getTrack(temp);

            cv::Scalar centroid = temp.at(0).getCentroidScalar();

            int x0 = (int) centroid[0];
            int y0 = (int) centroid[1];
            cv::Point p(y0,x0);
            color = cv::Scalar(rand() % 255, rand() % 255, rand() % 255);

            cv::circle(image, p, 1, color, -1, 8, 0);

            for(regit = temp.begin()+1; regit != temp.end(); regit++){
                
                cv::Scalar c = regit->getCentroidScalar();
                int x = (int) c[0];
                int y = (int) c[1];

                if(!isnan((double)(c[0])) && !isnan((double) c[1])){
                    cv::Point q(y,x);

                    cv::circle(image, q, 1, color, -1, 8, 0);
                    cv::line(image, p, q, color, 1);

                    p = q;
                    x0 = x;
                    y0 = y;
                }
            }
        }

        cv::imwrite(outfile, image);

    }

    /**
     * This function saves the complete set of cell tracks to one file.
     * For each track the following parameters are saved:
     * time point t
     * x-coordinate x
     * y-coordinate y
     * ID
     * boolean value for interaction (0/1)
     * area A
     *
     * @param file output filename
     * @param tracks vector of cell tracks
     */
    void printDataSet(const std::string file, std::vector<CellTrack> &tracks){

        std::vector<CellTrack>::iterator trackit;
	
        std::ofstream outfile(file);
	
        // build header
        if(outfile.is_open()){
            // outfile << "t\tx\ty\tID\tinteraction\tA\tflat\tMin\tMaj" << std::endl;
            // TODO: das hier mit stringstreams sauber machen
            outfile << "t\ty\tx\tID\tinteraction\tA\tflat\tMin\tMaj" << std::endl; // changed by Philipp: swap x and y
        }

        // fill data in file
        for(trackit = tracks.begin(); trackit != tracks.end(); trackit++){
            
            int t = trackit->getStartTime(); //get initial time
            std::vector<Region> regions;
            std::vector<Region>::iterator regit;
            
            trackit->getTrack(regions);
            int ID = trackit->getID();

            for(regit = regions.begin(); regit != regions.end(); regit++){
                cv::Scalar c = regit->getCentroidScalar();
                int interaction = trackit->getInteraction(t);

                if(outfile.is_open()){
                    outfile << t++ << "\t"; //increase time for the next step
                    outfile << c[1] << "\t";
                    outfile << c[0] << "\t";
                    outfile << ID << "\t";
                    outfile << interaction << "\t";
                    outfile << regit->getArea() << "\t";
                    
                    if(regit->getFlat())
                        outfile << "1" << "\t";
                    else
                        outfile << "0" << "\t";
                                        
                    outfile << regit->getMin()<< "\t";
                    outfile << regit->getMaj();
                    outfile << std::endl;
                }
            }
        }

        outfile.close();
    }

    /**
    * This function saves the complete set of regions within all cell tracks.
    *
    * @param path output filename
    * @param tracks of cell tracks
    */
    void printDataSetRegions(const std::string path, std::vector<CellTrack> &tracks){

        std::vector<CellTrack>::iterator trackit;

        std::stringstream ss;
        ss << path << "/regions.txt";

        std::string file = ss.str();
        std::ofstream outfile;
        outfile.open(file);

        if(outfile.is_open()){
            // TODO: das hier mit stringstream machen
            outfile << "ID\tt\tx\ty" << std::endl;
        }

        for(trackit = tracks.begin(); trackit != tracks.end(); trackit++){

            // get initial time
            int t = trackit->getStartTime();
            std::vector<Region> regions;
            std::vector<Region>::iterator regit;

            trackit->getTrack(regions);
            int ID = trackit->getID();

            for(regit = regions.begin(); regit != regions.end(); regit++){

                std::vector<cv::Vec2i> pixels = regit->region_pixels;

                for(int i = 0; i < pixels.size(); i++){
                    if(outfile.is_open()){
                        outfile << ID << "\t";
                        outfile << t << "\t"; // increase time for the next step
                        outfile << pixels.at(i).val[0] << "\t";
                        outfile << pixels.at(i).val[1];
                        outfile << std::endl;
                    }
                }

                t++;

            }
        }

        outfile.close();

    }

    /**
    * This function saves the all cluster regions over all gives time frames
    *
    * @param path output filename
    * @param tracks all tracks where cluster are included
    */
    void printClusterRegions(const std::string path, std::vector<CellTrack> &tracks){

        std::stringstream ss;
        ss << path << "/regions_cluster.txt";

        std::string file = ss.str();
        std::ofstream outfile;
        outfile.open(file);

        /// write headline
        if(outfile.is_open()){
            std::stringstream ss_headline;
            ss_headline << "ID" << "\t" << "t" << "\t" << "x" << "\t" << "y";
            outfile << ss_headline.str() << std::endl;
        }

        std::vector<CellTrack>::iterator trackit;

        /// iterate over all time frames
        for(trackit = tracks.begin(); trackit != tracks.end(); trackit++){

            // get initial time
            int t = trackit->getStartTime();
            std::vector<Region> regions;
            std::vector<Region>::iterator regit;

            trackit->getTrack(regions);
            int ID = trackit->getID();

            for(regit = regions.begin(); regit != regions.end(); regit++){

                /// check for cluster and extract info
                if ( regit->getClass() == cell_class::immune::CLUSTER ) {

                    std::vector<cv::Vec2i> pixels = regit->region_pixels;

                    for(int i = 0; i < pixels.size(); i++){
                        if(outfile.is_open()){
                            std::stringstream ss_content;
                            ss_content << regit->getId() << "\t" << t << "\t" << pixels[i].val[0] << "\t" << pixels[i].val[1];

                            outfile << ss_content.str() << std::endl;
                        }
                    }

                }

                t++;
            }
        }

        outfile.close();

    }

    /**
    * Draw an image on the previous segmented regions for class <RegionP>
    *
    * @param region_over_time regions
    * @param img_size underlying image size
    */
    void createImageFromRegion(cv::Mat &dst, std::vector<RegionP> &region_over_time, const cv::Size img_size, const cv::Scalar color, std::optional<int> thickness){

        cv::Mat img_mask = cv::Mat::zeros(img_size, CV_8UC1);

        std::vector<std::vector<cv::Point>> contours;

        // create all regions per frames
        std::vector<RegionP>::iterator reg;
        // iterate about all regions
        #pragma omp for firstprivate(region_over_time)
        for(reg = region_over_time.begin(); reg != region_over_time.end(); reg++){

            std::vector<cv::Point> cnt;

            // iterate over all contour points parse the cv::Vec2i contour vector into cv::Point
            for (size_t p = 0; p < reg->region_pixels.size(); ++p) {

                try {
                    cv::Vec2i region_point = reg->region_pixels[p];

                    cv::Point cnt_point = static_cast<cv::Point>(region_point);

                    cnt.push_back( cnt_point );
                } catch (const std::exception& e) {
                    std::cout << e.what() << std::endl;
                }

            }

            contours.push_back( cnt );

        }

        // draw the contour of the region into the image
        for(size_t c = 0; c < contours.size(); ++c ) {
            cv::drawContours(img_mask, contours, c, cv::Scalar(255), cv::FILLED);
        }

        // transpose the mask because of i-j differs for cv::Vec2i and cv::Point format
        cv::transpose(img_mask, img_mask);

        // some objects contours are not completely filled, find again all contours and fill them up
        std::vector<std::vector<cv::Point>> contours_2nd;
        std::vector<cv::Vec4i> hierarchy;
        cv::findContours( img_mask, contours_2nd, hierarchy, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_NONE, cv::Point(0, 0) );

        for(size_t c = 0; c < contours_2nd.size(); ++c ) {
            cv::drawContours(dst, contours_2nd, c, color, thickness.value_or(2), cv::LINE_8, hierarchy);
        }

    }

    /**
   * Draw an image on the previous segmented regions for class <Region>
   *
   * @param region_over_time regions
   * @param img_size underlying image size
   */
    void createImageFromRegion(cv::Mat &dst, std::vector<Region> &region_over_time, const cv::Size img_size, const cv::Scalar color, std::optional<int> thickness){

        cv::Mat img_mask = cv::Mat::zeros(img_size, CV_8UC1);

        std::vector<std::vector<cv::Point>> contours;

        // create all regions per frames
        std::vector<Region>::iterator reg;
        // iterate about all regions
        for(reg = region_over_time.begin(); reg != region_over_time.end(); reg++){

            std::vector<cv::Point> cnt;

            // iterate over all contour points parse the cv::Vec2i contour vector into cv::Point
            for (size_t p = 0; p < reg->region_pixels.size(); ++p) {

                try {
                    cv::Vec2i region_point = reg->region_pixels[p];

                    cv::Point cnt_point = static_cast<cv::Point>(region_point);

                    cnt.push_back( cnt_point );
                } catch (const std::exception& e) {
                    std::cout << e.what() << std::endl;
                }

            }

            contours.push_back( cnt );

        }

        // draw the contour of the region into the image
        for(size_t c = 0; c < contours.size(); ++c ) {
            cv::drawContours(img_mask, contours, c, cv::Scalar(255), cv::FILLED);
        }

        // transpose the mask because of i-j differs for cv::Vec2i and cv::Point format
        cv::transpose(img_mask, img_mask);

        // some objects contours are not completely filled, find again all contours and fill them up
        std::vector<std::vector<cv::Point>> contours_2nd;
        std::vector<cv::Vec4i> hierarchy;
        cv::findContours( img_mask, contours_2nd, hierarchy, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_NONE, cv::Point(0, 0) );

        for(size_t c = 0; c < contours_2nd.size(); ++c ) {
            cv::drawContours(dst, contours_2nd, c, color, thickness.value_or(2), cv::LINE_8, hierarchy);
        }

    }


} // outputs