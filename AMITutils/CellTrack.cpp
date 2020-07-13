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

#include <fstream>
#include "CellTrack.h"


/**
 * constructor of a new CellTrack object
 * @param t0 starting time
 * @param id ID of the CellTrack object
 */
CellTrack::CellTrack(int t0, int id) {
	this->t0 = t0;
	this->ID = id;
}

CellTrack::~CellTrack() {}

/**
 * This function adds a region to the end of the current CellTrack
 * @param region Object of type Region that is added to the current CellTrack
 */
void CellTrack::add(Region &region){
	region.setId(this->ID);
	this->track.push_back(region);
}

/**
 * This functions adds a new Region to the cell track at time point t.
 *
 * Timepoints can be in the range of [t_start-1, t_end +1].
 * If the track has already a measurement at timepoint t, the region will be replaced.
 *
 * @param region
 * @param t timepoint at which the region should be added.
 *
 */
void CellTrack::add(Region &region, const int &t){
	// add region at beginning of track
	if(t == this->t0 - 1){
		region.setId(this->ID);
		track.insert(track.begin(), region);
		this->t0--;
	}
	else if((unsigned) t == this->t0 + this->track.size()){
			track.insert(track.end(), region);
	}
	else if( t < t0){
		for(int i = t0 - 1 ; i > t; i--){
			Region r;
			r.setCentroid(NAN, NAN);
			track.insert(track.begin(),r);
		}
		track.insert(track.begin(), region);
	}
	else if ((unsigned) t  > this->t0 + this->track.size()){
		for(int i =  this->t0 + this->track.size(); i < t; i++){
			Region r;
			r.setCentroid(NAN, NAN);
			track.insert(track.end(),r);
		}
		track.insert(track.end(), region);
	}
	// t is neither t_start-1 nor t_end+1! --> cannot add region
	else{ 
		this->change(region, t);
	}
}

/**
 * This function adds a vector of regions to the end of the current CellTrack
 * @param regions vector of regions to be added to the object
 */
void CellTrack::add(std::vector<Region> &regions){
	std::vector<Region>::iterator it;
	int id = this->ID;
	
    for(it = regions.begin(); it != regions.end(); it++){
		it->setId(id);
		this->add((*it));
	}
}

/**
 * This function adds another CellTrack to the current CellTrack.
 *
 * Starting and ending times are set according to the minimum and maximum timepoints of both tracks, respectively.
 * Missing timepoints are set to missing values (centroid points of [NA,NA]).
 * Region pixels are merged and contour pixels are recalculated.
 *
 *
 * @param track1 CellTrack that is added to the object
 */
void CellTrack::add(CellTrack &track1){
	std::vector<Region> tmp;
	track1.getTrack(tmp);

	int tstart = track1.getStartTime();
	if(tstart == this->getEndTime() + 1){
		this->add(tmp);
	}

	else if(tstart > this->getEndTime() + 1){
		// add nan regions until tstart is reached
		int dt = tstart - this->getEndTime() - 1;
		for(size_t i = 0; i < dt; i++){
			Region r;
			r.setCentroid(NAN, NAN);
			this->add(r);
		}

		this->add(tmp);
	}
	// overlap in time
	else if(tstart <= this->getEndTime()){

		if(this->getStartTime() > tstart){
			// add regions to beginning of the track;
			for(size_t t =  this->getStartTime()-1; t >= tstart; t--){
				Region r;
				track1.getRegionAt(t,r);
				this->add(r, t);
			}
		}

		for(size_t i = this->getStartTime(); i <= this->getEndTime(); i++){
			Region r;
			track1.getRegionAt(i, r);
			track.at(i-this->t0).addRegionPixels(r.region_pixels);
			track.at(i-this->t0).computeContour();
			track.at(i-this->t0).computeCentroid();
		}

		for(size_t i = this->getEndTime()+1; i <= track1.getEndTime(); i++){
			Region r;
			track1.getRegionAt(i, r);
			this->add(r);
		}
	}
}

/**
 * This function alters the Region at time point t.
 *
 * @param region
 * @param t
 */
void CellTrack::change(const Region &region, const int &t){
	if(t >= this->t0 && t <= this->getEndTime()){
		this->track.at(t-this->t0) = region;
	}
}

/**
 * This function checks if the measurement at timepoint t is valid (no missing measurement)
 *
 * @param t timepoint
 *
 * @return true if track has no missing measurement at timepoint t
 */
bool CellTrack::hasValidMeasurement(const int &t){
	if(t < this->t0 || t > this->getEndTime())
		return false;
	
	cv::Point center =  this->track.at(t-this->t0).getCentroidPoint();
	if(center.x >= 0 && center.y >= 0)
		return true;
	else
        return false;
}

/**
 * This function checks if the CellTrack exists at specific timepoint.
 *
 * @param t timepoint
 *
 * @return true if the cell track includes timepoint t
 */
bool CellTrack::existsAt(const int &t){
	if( (t >= this->t0)  && (t < (this->t0 + this->getLength()))){
		return true;
	}
	return false; //else
}

/**
 * This function returns the ID of the current cell track.
 *
 * @return ID
 */
int CellTrack::getID(){
	return this->ID;
}

/**
 * This function returns the length of the current CellTrack
 *
 * @return track length
 */
int CellTrack::getLength(){
	return this->track.size();
}

/**
 * @return starting time
 */
int CellTrack::getStartTime(){
	return this->t0;
}

/**
 * This function returns the ending time of the CellTrack.
 *
 * @return ending time
 */
const int CellTrack::getEndTime(){
	return (this->t0 + this->getLength() - 1);
}

/**
 * This function gives access to the region at a specific timepoint.
 * This is a call-by-reference, changes made to the parameter affect the passed argument!
 *
 * @param t timepoint
 * @param region
 */
void CellTrack::getRegionAt(const int &t, Region &region){
	if( this->existsAt(t)){
		int index = t - this->t0;
		region = this->track.at(index);
	}
}

/**
 * This function gives access to the region at a specific timepoint.
 * This is a call-by-reference, changes made to the parameter affect the passed argument!
 *
 * @param t timepoint
 * @return region by address
 */
Region* CellTrack::getRegionAt(const int &t){
	if( this->existsAt(t)){
		int index = t - this->t0;
		return &this->track.at(index);
	}
	else{
		return NULL;
	}
}

/**
 * This function computes the average area of the current CellTrack.
 *
 * @return average area
 */
int CellTrack::getAverageArea(){
	int area = 0;

	for(size_t i = 0; i < this->track.size(); i++){
		if(this->hasValidMeasurement(i+this->t0)){
			area += this->track.at(i).getArea();
		}
	}

	area /= this->getLength();
	
	return area;
}

/**
 * This function computes the maximum speed of the current CellTrack.
 *
 * @return maximum speed value
 */
double CellTrack::getMaxSpeed(){
	double maxspeed = 0.0;
	std::vector<Region>::iterator it, it2;

	if(this->track.size() < 2)
		return 0.0;
	
	for(it = this->track.begin(), it2 = this->track.begin()+1; it != this->track.end()-1, it2 != this->track.end(); it++, it2++){
		cv::Scalar c1 = it->getCentroidScalar();
		cv::Scalar c2 = it2->getCentroidScalar();

		double x1 = c1[0];
		double x2 = c2[0];
		double y1 = c1[1];
		double y2 = c2[1];

		double speed = sqrt( pow((x2-x1), 2.) + pow((y2-y1), 2.) );
		if(speed > maxspeed){
			maxspeed = speed;
		}
	}

	return maxspeed;
}

/**
 * This function gives access to the first region.
 * This is a call-by-reference, changes made to the parameter affect the passed argument!
 *
 * @param region
 */
void CellTrack::getFirstRegion(Region &region){
	if(track.size() > 0)
		region = this->track.at(0);	
	else
		std::cout << "Error! Track size zero! "<< std::endl;	
}

/**
 * This function gives access to the last region.
 * This is a call-by-reference, changes made to the parameter affect the passed argument!
 *
 * @param region
 */
void CellTrack::getLastRegion(Region &region){
	if(track.size() > 0)
		region = this->track.at(track.size()-1);	
	else
		std::cout << "Error! Track size zero! "<< std::endl;	
}

/**
 * This function gains access to the track.
 *
 * This is a call-by-reference, changes made to the parameter affect the passed argument!
 *
 * @param track (vector of type Region)
 */
void CellTrack::getTrack(std::vector<Region> &track){
	track = this->track;
}

/**
 * This function returns the length of missing measurements of the current cell track
 *
 * @return number of missing measurements
 */
int CellTrack::getNumberOfMissingMeasurements(){
	int n = 0;

	std::vector<Region>::iterator it;
	for(it = this->track.begin(); it != this->track.end(); it++){
		cv::Scalar c = it->getCentroidScalar();
		
		if(isnan(c[0]) || isnan(c[1])){
			n++;
		}
	}
	
	return n;
}

/**
 * This function computes the number of interactions with other cells.
 *
 * @return number of interactions
 */
int CellTrack::getNumberOfInteractions(){
	std::vector<Region>::iterator regit;
	
	int n = 0;
	for(regit =  this->track.begin(); regit!= this->track.end(); regit++){
		if(regit->isInteracting()){
			n++;
		}
	}
	
	return n;
}

bool CellTrack::getInteraction(const int &t){
	return this->track.at(t-t0).isInteracting();
}

/**
 * This function computes the minimum enclosing ellipse of the region at timepoint t.
 *
 * @param t timepoint
 * @param size ellipse size
 */
void CellTrack::getPrincipalAxesMinEnclosingEllipse(int t, cv::Size2f &size){
	this->track.at(t-this->t0).getPrincipalAxesMinEnclosingEllipse(size);
}

/**
 * This function sets an interaction event for time point t and cell ID.
 *
 * @param t time point
 * @param id fungal cell ID
 */
void CellTrack::setInteraction(int t, int id){
	this->track.at(t-this->t0).addInteractionId(id);
}

/**
 * This function sets the cell track ID.
 *
 * @param id ID
 */
void CellTrack::setId(int id){
	this->ID = id;

	for(std::vector<Region>::iterator regit = this->track.begin(); regit != this->track.end(); regit++){
		regit->setId(id);
	}
}

/**
 * This function sets the start time of the cell track.
 *
 * @param t start time
 */
void CellTrack::setStartTime(int t){
	this->t0 = t;
}

void CellTrack::clearOverlaps(){
	this->overlapregionsb1.clear();
	this->overlapregionsb2.clear();
	this->overlapregionsf1.clear();
	this->overlapregionsf2.clear();

	this->overlapregionsb1_fc.clear();
	this->overlapregionsb2_fc.clear();
	this->overlapregionsf1_fc.clear();
	this->overlapregionsf2_fc.clear();
}

/**
 * This function changes interaction IDs (i.e. IDs of interacting cells)
 * @param id old ID
 * @param id_new new ID
 */
void CellTrack::changeInteractionID(int id, int id_new){
	
	std::vector<Region>::iterator it = this->track.begin();
	for(; it != this->track.end(); it++){
		it->changeInteractionID(id, id_new);
	}
}

/**
 * This function shortens a track until timepoint t (excluding t).
 * Therefore all track measurements after t are removed.
 *
 * @param t timepoint
 */
void CellTrack::removeRegions(int t){
	this->track.resize(t-t0+1);
}

/**
 * This function saves the CellTrack to a file.
 * The output includes time point, x- and y-coordinates, ID, area, boolean value for interaction and interaction IDs
 *
 *
 * @param path file path
 */
void CellTrack::printToFile(const std::string path){
	std::stringstream ss;
	ss << path << "/" << this->getID() << ".txt";
	
	std::string file = ss.str();	
	std::ofstream outfile(file);
	
	// build headline
	if(outfile.is_open()){
		outfile << "t\tx\ty\tID\tA\tI\tI_ids\tn_regions\tflat\tMin\tMaj" << std::endl;
		// outfile << "t\ty\tx\tID\tA\tI\tI_ids\tn_regions\tflat\tMin\tMaj" << std::endl;  // changed by Philipp: swap x and y
	}

	int t = this->getStartTime();
	std::vector<Region>::iterator it;

	// enter data
	for(it = this->track.begin(); it != this->track.end(); it++){
		
		cv::Scalar centroid = it->getCentroidScalar();
		double x = centroid[1];
		double y = centroid[0];

		if(!std::isnan(x) && !std::isnan(y)){
			outfile << t << "\t" << x << "\t" << y << "\t" << this->ID;
			outfile << "\t" << it->getArea();

			if(it->isInteracting())
				outfile << "\t" << 1;			
			else
				outfile << "\t" << 0;			

			// interaction ids
			std::vector<int> ids;
			it->getInteractionIds(ids);
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

			outfile << "\t";
			outfile << it->n_regions;

			outfile << "\t";
			
			if(it->getFlat())
				outfile << "1";
			else
				outfile << "0";
			
			if(it->getArea() > 0){
				outfile << "\t";
				outfile << it->getMin();
				outfile << "\t";
				outfile << it->getMaj();
			}else{
				outfile << "\t0/t0";
			}
		}
		else{
			outfile << t << "\t" << "NA\tNA\t" << this->ID<< "\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
		}

		outfile << std::endl;
		t++;
	}
	outfile.close();
}
