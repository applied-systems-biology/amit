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

#include "CellTrackP.h"
#include <fstream>

CellTrackP::CellTrackP(int t0, int id) : CellTrack(t0, id) {}

CellTrackP::~CellTrackP() {}

/**
 * This function adds a region to the end of the current CellTrack
 * @param region Object of type Region that is added to the current CellTrack
 */
void CellTrackP::add(RegionP &region){
	region.setId(this->ID);
	this->track.push_back(region);
}

/**
 * This function adds a vector of regions to the end of the current CellTrack
 * @param regions vector of regions to be added to the object
 */
void CellTrackP::add(std::vector<RegionP> &regions){
	std::vector<RegionP>::iterator it;
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
void CellTrackP::add(CellTrackP &track1){
	std::vector<RegionP> tmp;
	track1.getTrack(tmp);

	int tstart = track1.getStartTime();
	if(tstart == this->getEndTime() + 1){
		this->add(tmp);
	}

	else if(tstart > this->getEndTime() + 1){
		//add nan regions until tstart is reached
		int dt = tstart - this->getEndTime() - 1;
		for(int i = 0; i < dt; i++){
			RegionP r;
			r.setCentroid(NAN, NAN);
			this->add(r);
		}
		this->add(tmp);

	}

	//overlap in time
	else if(tstart <= this->getEndTime()){

		if(this->getStartTime() > tstart){
			//add regions to beginning of the track;
			for(int t =  this->getStartTime()-1; t >= tstart; t--){
				RegionP r;
				track1.getRegionAt(t,r);
				this->add(r, t);
			}
		}

		for(int i = this->getStartTime(); i <= this->getEndTime(); i++){
			RegionP r;
			track1.getRegionAt(i, r);
			track.at(i-this->t0).addRegionPixels(r.region_pixels);
			track.at(i-this->t0).computeContour();
			track.at(i-this->t0).computeCentroid();
		}

		for(int i = this->getEndTime()+1; i <= track1.getEndTime(); i++){
			RegionP r;
			track1.getRegionAt(i, r);
			this->add(r);
		}
	}

//	//remove phagocytoses for second cell track
//	int t = max(track1.getStartTime(), this->getStartTime());
//	cout << t << endl;
//	for(; t <= this->getEndTime(); t++){
//		cout << t << endl;
//		if(!this->hasValidMeasurement(t)){
//			continue;
//		}
//		cout << track.at(t-t0).getCentroidPoint() << endl;
//		if(this->track.at(t-t0).getNumberOfPhagocytoses() > 0){
//			cout << "\tvalid" << endl;
//
//			this->track.at(t-t0).removePhIDs();
//			cout << "\tremoved" << endl;
//		}
//	}
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
 *
 */
void CellTrackP::add(RegionP &region, int t){
	//add region at beginning of track
	if(t == this->t0 -1){
		region.setId(this->ID);
		track.insert(track.begin(), region);
		this->t0--;
	}
	else if((unsigned) t == this->t0 + this->track.size()){
			track.insert(track.end(), region);
	}
	else if( t < t0){
		for(int i = t0-1 ; i > t; i--){
			RegionP r;
			r.setCentroid(NAN, NAN);
			track.insert(track.begin(),r);
		}
		track.insert(track.begin(), region);
	}
	else if ((unsigned) t  > this->t0 + this->track.size()){
		for(int i =  this->t0 + this->track.size(); i < t; i++){
			RegionP r;
			r.setCentroid(NAN, NAN);
			track.insert(track.end(),r);
		}
		track.insert(track.end(), region);
	}
	else{ //t is neither t_start-1 nor t_end+1! --> cannot add region
		this->change(region, t);
//		cout << "Error: Add region to track " << this->ID << " at time " << t << " failed!" << endl;
	}
}

void CellTrackP::addFungiId(int t, int id){
	this->track.at(t-this->t0).addFungiId(id);
}

void CellTrackP::addPhagocytosis(int tstart, int id){
	int pos =  this->track.at(tstart-this->t0).getNumberOfPhagocytoses();
//	cout << "pos " << pos << endl;
	for(; tstart <= this->getEndTime(); tstart++){
		if(this->hasValidMeasurement(tstart)){
//			cout << tstart << endl;
			this->track.at(tstart - this->t0).addPhagocytosisIDPos(id, pos);

		}

	}
}

void CellTrackP::addPhagocytosis(int tstart, int id, int pos){
//	int pos =  this->track.at(tstart-this->t0).getNumberOfPhagocytoses();
//	cout << "pos " << pos << endl;
	for(; tstart <= this->getEndTime(); tstart++){
		if(this->hasValidMeasurement(tstart)){
//			cout << tstart << endl;
			this->track.at(tstart - this->t0).addPhagocytosisIDPos(id, pos);

		}

	}
}

/**
 * This function alters the Region at time point t.
 *
 * @param region
 * @param t
 */
void CellTrackP::change(RegionP & region, int t){
	if(t >= this->t0 && t <= this->getEndTime()){
		this->track.at(t-this->t0) = region;
	}
}

void CellTrackP::changeFID(int id, int id_new){
	for(unsigned int t = 0; t <  this->track.size(); t++){
		this->track.at(t).changeFID(id, id_new);
		this->track.at(t).changePhagocytosisID(id, id_new);
	}
}

/**
 * This function deletes the last track region.
 */
void CellTrackP::deleteLastRegion(){
	this->track.pop_back();
}

void CellTrackP::clearOverlaps(){
	this->overlapregionsb1.clear();
	this->overlapregionsb2.clear();
	this->overlapregionsf1.clear();
	this->overlapregionsf2.clear();
}

/**
 * This function deletes the region at a specific timepoint.
 * If the timepoint is starting or ending time, the track length is reduced by 1.
 * Else, the centroid replaced by missing values.
 *
 *@param t timepoint
 */
void CellTrackP::deleteRegion(int t){
	if(t > this->t0 && t < this->getEndTime()){
		RegionP r;
		r.setCentroid(NAN, NAN);
		this->track.at(t - this->t0) = r;
	}
	else if(t == this->t0){
		this->deleteFirstRegion();
	}
	else if(t == this-> getEndTime()){
		this->deleteLastRegion();
	}
}

/**
 * This function deletes the first track region.
 */
void CellTrackP::deleteFirstRegion(){
	this->track.erase(this->track.begin());
	this->t0++;
}


bool CellTrackP::existsAt(int t){
	if( (t >= this->t0)  && (t < (this->t0 + this->getLength()))){
		return true;
	}
	return false; //else
}

/**
 * This function computes the average area of the current CellTrack.
 *
 * @return average area
 */
int CellTrackP::getAverageArea(){
	int area = 0;

	for(unsigned int i = 0; i < this->track.size(); i++){
		area += (double) this->track.at(i).getArea();
	}
	area /= this->track.size();
	return area;
}

/**
 * This function returns the ending time of the CellTrack.
 *
 * @return ending time
 */
const int CellTrackP::getEndTime(){
	return (this->t0 + this->getLength() - 1);
}

void CellTrackP::getFIDs(int t, std::vector<int> &fids){
	this->track.at(t-this->t0).getFungiIds(fids);
}


/**
 * This function gives access to the first region.
 * This is a call-by-reference, changes made to the parameter affect the passed argument!
 *
 * @param region
 */
void CellTrackP::getFirstRegion(RegionP &region){
	if(track.size() > 0){
		region = this->track.at(0);
	}
	else{
		std::cout << "Error! Track size zero! "<< std::endl;
	}
}

/**
 * This function returns the timepoint of first phagocytosis
 * If no phagocytosis occured, -1 is returned
 * @return
 */
int CellTrackP::getFirstPhagocytosisTimepoint(){
	int t = this->t0;

	while(t <= this->getEndTime() && this->getNumberOfPhagocytoses(t) < 1){
		t++;
	}

	if(t > this->getEndTime()){
		return -1;
	}

	else{
		return t;
	}
}

/**
 * This function checks if the measurement at timepoint t is valid (no missing measurement)
 *
 * @param t timepoint
 *
 * @return true if track has no missing measurement at timepoint t
 */
bool CellTrackP::hasValidMeasurement(int t){

	if(t < this->t0 || t > this->getEndTime()){
		return false;
	}
	cv::Point center =  this->track.at(t-this->t0).getCentroidPoint();
	if(center.x >= 0 && center.y >= 0){
		return true;
	}
	//else
	return false;
}


int CellTrackP::getIState(int t){
	return this->track.at(t-this->t0).getIState();
}

void CellTrackP::getIStates(std::vector<int> &states){

	states.resize(0);
	std::vector<RegionP>::iterator regit;
	for(regit = this->track.begin(); regit != this->track.end(); regit++){
		states.push_back(regit->getIState());
	}

}

/**
 * This function gives access to the last region.
 * This is a call-by-reference, changes made to the parameter affect the passed argument!
 *
 * @param region
 */
void CellTrackP::getLastRegion(RegionP &region){
	if(track.size() > 0){
		region = this->track.at(track.size()-1);
	}
	else{
		std::cout << "Error! Track size zero! "<< std::endl;
	}
}

/**
 * This function returns the length of the current CellTrack
 *
 * @return track length
 */
int CellTrackP::getLength(){
	return this->track.size();
}

/**
 * This function computes the maximum speed of the current CellTrack.
 *
 * @return maximum speed value
 */
double CellTrackP::getMaxSpeed(){
	double maxspeed = 0.0;
	std::vector<RegionP>::iterator it, it2;

	if(this->track.size() < 2){
		return 0.0;
	}
//	cout << this->track.size() << endl;

	for(it = this->track.begin(), it2 = this->track.begin()+1; it != this->track.end()-1, it2 != this->track.end(); it++, it2++){
		cv::Scalar c1 = it->getCentroidScalar();
		cv::Scalar c2 = it2->getCentroidScalar();
//		cout << c1 << "," << c2 << endl;

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

int CellTrackP::getNumberOfFungi(int t){
	std::vector<int> ids;
	this->track.at(t-this->t0).getFungiIds(ids);
	return ids.size();
}

int CellTrackP::getNumberOfIState(int i){
	std::vector<RegionP>::iterator regit;
	int n = 0;
	for(regit =  this->track.begin(); regit!= this->track.end(); regit++){
		if(regit->getIState() == i){
			n++;
		}
	}
	return n;
}

/**
 * This function returns the length of missing measurements of the current cell track
 *
 * @return number of missing measurements
 */
int CellTrackP::getNumberOfMissingMeasurements(){

	int n = 0;
	std::vector<RegionP>::iterator it;
	for(it = this->track.begin(); it != this->track.end(); it++){
		cv::Scalar c = it->getCentroidScalar();
		if(std::isnan(c[0]) || std::isnan(c[1])){
			n++;
		}
	}
	return n;
}

int CellTrackP::getNumberOfPhagocytoses(int t){
	return this->track.at(t-this->t0).getNumberOfPhagocytoses();
}

int CellTrackP::getNumberOfPhagocytozedFungi(int t){
	return this->track.at(t-this->t0).getNumberOfPhagocytoses();
}

void CellTrackP::getPhagocytosisIds(int t, std::vector<int> &ids){
	this->track.at(t-this->t0).getPhagocytosisIds(ids);
}

/**
 * This function gives access to the region at a specific timepoint.
 * This is a call-by-reference, changes made to the parameter affect the passed argument!
 *
 * @param t timepoint
 * @param region
 */
RegionP * CellTrackP::getRegionAt(int t){
	if( this->existsAt(t)){
		int index = t - this->t0;
		return &(this->track.at(index));
	}
	else return NULL;
}

/**
 * This function gives access to the region at a specific timepoint.
 * This is a call-by-reference, changes made to the parameter affect the passed argument!
 *
 * @param t timepoint
 * @param region
 */
void CellTrackP::getRegionAt(int t, RegionP & region){
	if( this->existsAt(t)){
		int index = t - this->t0;
		region = this->track.at(index);
	}
}

void CellTrackP::getTimepointsOfPhagocytoses(std::vector<int> & timepoints){
	timepoints.clear();
	int n_ph = 0;

	int t = this->t0;
	for(std::vector<RegionP>::iterator regit = this->track.begin(); regit != this->track.end(); regit++, t++){
		if(regit->getNumberOfPhagocytoses() > n_ph){
			timepoints.push_back(t);
			n_ph++;
		}
	}
}

/**
 * This function gains access to the track.
 *
 * This is a call-by-reference, changes made to the parameter affect the passed argument!
 *
 * @param track (vector of type Region)
 */
void CellTrackP::getTrack(std::vector<RegionP> &track){
	track = this->track;
}

bool CellTrackP::hasFungiId(int id){
	for(int t = 0; t <= this->getEndTime()-this->t0; t++){
		if(this->track.at(t).hasFungiID(id)){
			return true;
		}
	}
	//else
	return false;
}


bool CellTrackP::hasFungiId(int id, int t){
	if(this->track.at(t-this->t0).hasFungiID(id)){
		return true;
	}
	//else
	return false;
}

bool CellTrackP::hasPhagocytozedFungus(int id){
	return this->track.back().hasPhagocytozedFungus(id);
}

/**
 *
 * Determines whether the cell is interaction with another cell at a specific timepoint.
 *
 * @param t time point
 *
 * @return true, if CellTrack is interacting at time point t
 */
bool CellTrackP::isInteracting(int t){
//	cout << "t " << t << endl;
//	cout << this->track.size() << endl;
//	cout << "t-t0 " << t - this->t0 << endl;
	RegionP p = this->track.at(t-this->t0);
//	cout << p.getCentroidPoint() << endl;
	std::vector<int> ids;
	p.getInteractionIds(ids);
//	cout << ids.size() << endl;

	return this->track.at(t-this->t0).isInteracting();
}

bool CellTrackP::isTouchingMoreFungiThan(int t, int fid){
	RegionP r = track.at(t-this->t0);
	std::vector<int> interacting_fungi;

	r.getFungiIds(interacting_fungi);
	if(interacting_fungi.size() > 1){
		return true;
	}
	//else
	return false;
}

/**
 * This function saves the CellTrack to a file.
 * The output includes time point, x- and y-coordinates, ID, area, boolean value for interaction and interaction IDs
 *
 *
 * @param path file path
 */
void CellTrackP::printToFile(std::string path){
	
	std::stringstream ss;
	ss << path << this->getID() << ".txt";
	std::string file = ss.str(); 

	std::ofstream outfile;
	outfile.open(file);

	if(outfile.is_open()){
		outfile << "t\tx\ty\tID\tFR\tA\tS\tnPh\tI\tI_ids\tF_ids\tDF_ids\tPh_ids\tState_probabilities" << std::endl;
	}

	int t = this->getStartTime();
	std::vector<RegionP>::iterator it;

	for(it = this->track.begin(); it != this->track.end(); it++){
		cv::Scalar centroid = it->getCentroidScalar();
		double x = centroid[1];
		double y = centroid[0];

		if(!std::isnan(x) && !std::isnan(y)){
			outfile << t << "\t" << x << "\t" << y << "\t" << this->ID << "\t" << it->getFungalRatioGreen();
			outfile << "\t" << it->getArea();

			outfile << "\t" << it->getIState();
			outfile << "\t" << it->getNumberOfPhagocytoses();
			if(it->isInteracting()){
				outfile << "\t" << 1;
			}
			else{
				outfile << "\t" << 0;
			}

			/// interaction ids
			std::vector<int> ids;
			it->getInteractionIds(ids);
			if(ids.size() == 0){
				outfile << "\t" << 0;
			}
			else{
				outfile << "\t";
				for(unsigned int i = 0; i < ids.size()-1; i++){
					outfile << ids.at(i) << ",";
				}
				outfile << ids.at(ids.size()-1);
			}

			/// fungi ids
			it->getFungiIds(ids);
			outfile << "\t";
			if(ids.size() == 0){
				outfile << 0;
			}
			else{
				for(unsigned int i = 0; i < ids.size()-1; i++){
					outfile << ids.at(i) << ",";
				}
				outfile << ids.at(ids.size()-1);
			}

			/// dead fungi ids
			it->getDeadFungiIds(ids);
			outfile << "\t";
			if(ids.size() == 0){
				outfile << 0;
			}
			else{
				for(unsigned int i = 0; i < ids.size()-1; i++){
					outfile << ids.at(i) << ",";
				}
				outfile << ids.at(ids.size()-1);
			}

			//phagocytosis ids
			it->getPhagocytosisIds(ids);
			outfile << "\t";
			if(ids.size() == 0){
				outfile << 0;
			}
			else{
				for(unsigned int i = 0; i < ids.size()-1; i++){
					outfile << ids.at(i) << ",";
				}
				outfile << ids.at(ids.size()-1);
			}

		}
		else{
			outfile << t << "\t" << "nan\tnan\t" << this->ID<< "\t0\t0\t0\t0\t0\t0\t0\t0\t0";
		}
		outfile << std::endl;
		t++;
	}
	outfile.close();
}

void CellTrackP::removeFungiIds(int t){
	track.at(t-this->t0).removeFungiIds();
}

void CellTrackP::removeFungiId(int t, int id){
	this->track.at(t-this->t0).removeFungiId(id);
}

void CellTrackP::removeLastPhID(int t){
//	cout << "remove at t " << t << endl;
	this->track.at(t-this->t0).removeLastPhID();
}

/**
 * This function removes missing values at beginning and end of the track.
 * The track is shortened accordingly.
 */
void CellTrackP::removeNaNRegions(){
//	cout << this->track.size() + t0 << endl;

	if(track.size() > 0){

		RegionP r;
		this->getLastRegion(r);
		cv::Point center = r.getCentroidPoint();

		while(this->getLength() > 0){
			center = r.getCentroidPoint();
			if((center.x <0 || center.y < 0)){
				this->deleteLastRegion();
				this->getLastRegion(r);
			}
			else{
				break;
			}
		}

		this->getFirstRegion(r);
		while(this->getLength() > 0){

			center = r.getCentroidPoint();
			if((center.x <0 || center.y < 0)){
				this->deleteFirstRegion();
				this->getFirstRegion(r);
			}
			else{
				break;
			}
		}

	}
}

/**
 * This function shortens a track until timepoint t (excluding t).
 * Therefore all track measurements after t are removed.
 *
 * @param t timepoint
 */
void CellTrackP::removeRegions(int t){
	this->track.resize(t-t0+1);
}

/**
 * This function sets an interaction event for time point t and cell ID.
 *
 * @param t time point
 * @param id fungal cell ID
 */
void CellTrackP::setInteraction(int t, int id){
//	this->track.at(t-this->t0).setInteraction(true);
	this->track.at(t-this->t0).addInteractionId(id);
}

/**
 * This function removes a region at a specific timepoint and replaces it by missing values.
 * If the time point is the ending or starting time, then the track will be shortened.
 *
 * @param t time point
 */
void CellTrackP::removeRegionAt(int t){
	if(t >= this->t0 && t <= this->getEndTime()){
		RegionP r;
		r.setCentroid(NAN, NAN);
		this->track.at(t - this->t0) = r;

		this->removeNaNRegions();
	}

}

/**
 *
 */
void CellTrackP::setIState(int t, int state){
	this->track.at(t-this->t0).setIState(state);
}
