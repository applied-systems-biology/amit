/*
 * RegionP.cpp
 *
 *  Created on: 04.08.2015
 *      Author: susanne
 */

#include "RegionP.h"
#include <opencv2/highgui.hpp>
#include <algorithm>
#include "Calcs.h"


RegionP::RegionP() : Region() {
	this->I_state = 0;
}

RegionP::RegionP(Region &r) : Region(r) {
	this->I_state = 0;
}

RegionP::RegionP(Region &r, bool fc) : Region(r) {
	this->I_state = 0;
	this->FC = fc;
}

RegionP::~RegionP() { }



/**
 * This function includes the ID of an interacting fungal region.
 * IDs cannot be added more than once.
 * The list of IDs is sorted in numerically increasing mode.
 */
void RegionP::addFungiId(int id){

	if(std::find(this->fungiIDs.begin(), this->fungiIDs.end(), id) == this->fungiIDs.end()){
		this->fungiIDs.push_back(id);
		sort(this->fungiIDs.begin(), this->fungiIDs.end());
	}
}

/**
 * This function includes the ID of a dead fungal region.
 * IDs cannot be added more than once.
 * The list of IDs is sorted in numerically increasing mode.
 */
void RegionP::addDeadFungiId(int id){

	if(std::find(this->deadFungiIDs.begin(), this->deadFungiIDs.end(), id) == this->deadFungiIDs.end()){
		this->deadFungiIDs.push_back(id);
		sort(this->deadFungiIDs.begin(), this->deadFungiIDs.end());
	}
}

/**
 * This function adds a single pixel to the set of contour pixels.
 */
void RegionP::addGreenPx(const cv::Vec2i &p){
	if(std::find(this->green_region_pixels.begin(), this->green_region_pixels.end(), p) == this->green_region_pixels.end()){
		this->green_region_pixels.push_back(p);
	}
}

void RegionP::addGreenPxs(std::vector<cv::Vec2i> & px){

	std::vector<cv::Vec2i>::iterator it;
	for(it = px.begin(); it != px.end(); it++){
		this->addGreenPx(*it);
	}
}

/**
 * This function adds the ID of a phagocytozed fungal cell to the region.
 * The vector is not sorted to save the chronological order of the phagocytosis events.
 */
void RegionP::addPhagocytosisID(int id){
	if(std::find(this->phagocytosisIDs.begin(), this->phagocytosisIDs.end(), id ) == this->phagocytosisIDs.end()){
		this->phagocytosisIDs.push_back(id);
	}
}

void RegionP::addPhagocytosisIDPos(int id, int pos){

	if(std::find(this->phagocytosisIDs.begin(), this->phagocytosisIDs.end(), id ) == this->phagocytosisIDs.end()){
//		for(unsigned int i = 0; i < this->phagocytosisIDs.size(); i++){
//			cout << this->phagocytosisIDs.at(i) << ", ";
//		}
//		cout << endl;
		this->phagocytosisIDs.insert( this->phagocytosisIDs.begin() + pos, id);
//		for(unsigned int i = 0; i < this->phagocytosisIDs.size(); i++){
//			cout << this->phagocytosisIDs.at(i) << ", ";
//		}
//		cout << endl;
	}
	else{
		std::vector<int>::iterator found = std::find(this->phagocytosisIDs.begin(), this->phagocytosisIDs.end(), id );
		std::vector<int>::iterator fpos = this->phagocytosisIDs.begin() + pos;
		if(found != fpos){
			this->phagocytosisIDs.erase(found);
			if(pos < (signed) this->phagocytosisIDs.size()){
				this->phagocytosisIDs.insert( this->phagocytosisIDs.begin() + pos, id);
			}
			else if(pos == (signed) this->phagocytosisIDs.size()){
				this->phagocytosisIDs.push_back(id);
			}
		}
	}
}


void RegionP::changeFID(int id_old, int id_new){
	if(std::find(this->fungiIDs.begin(), this->fungiIDs.end(), id_old) != this->fungiIDs.end()){

		std::vector<int>::iterator it = std::find(this->fungiIDs.begin(), this->fungiIDs.end(), id_new);

		//avoid double entries
		if(it != this->fungiIDs.end()){
			//remove id_old
			this->fungiIDs.erase(std::remove(this->fungiIDs.begin(), this->fungiIDs.end(), id_old));
		}
		else{
			std::replace(this->fungiIDs.begin(), this->fungiIDs.end(), id_old, id_new);
		}
	}
}

void RegionP::changePhagocytosisID(int id_old, int id_new){

	if(std::find(this->phagocytosisIDs.begin(), this->phagocytosisIDs.end(), id_old) != this->phagocytosisIDs.end()){

		std::vector<int>::iterator it = std::find(this->phagocytosisIDs.begin(), this->phagocytosisIDs.end(), id_new);

		//avoid double entries
		if(it != this->phagocytosisIDs.end()){
			//remove id_old
			this->phagocytosisIDs.erase(std::remove(this->phagocytosisIDs.begin(), this->phagocytosisIDs.end(), id_old));
		}
		else{
			std::replace(this->phagocytosisIDs.begin(), this->phagocytosisIDs.end(), id_old, id_new);
		}

	}
}

/**
 * This function computes the centroid of the region from all region pixels that are not present in the set of green region pixels.
 */
void RegionP::computeCentroidFromGrayPixels(){
	std::vector<cv::Vec2i>::iterator it;
	int n =  this->getArea() - this->green_region_pixels.size();

	if(n  > 0){
		double sumi = 0.0, sumj = 0.0;

		for(it = this->region_pixels.begin(); it != this->region_pixels.end(); it++){
			//only add region pixels that do not belong to the green pixel sets
			if(! Calcs::exists(*it, this->green_region_pixels)){
				sumi += (*it)[0];
				sumj += (*it)[1];
			}
		}
		sumi /= n;
		sumj /= n;

		this->centroid[0] = sumi;
		this->centroid[1] = sumj;
	}


}

/**
 * This function computes the center of mass from the region pixels that only belong to the set of green region pixels but not to the set of region pixels.
 */
void RegionP::computeCentroidFromGreenPixels(){
	std::vector<cv::Vec2i>::iterator it;
	int n =  this->green_region_pixels.size();

	if(n  > 0){
		double sumi = 0.0, sumj = 0.0;

		for(it = this->green_region_pixels.begin(); it != this->green_region_pixels.end(); it++){
			//only add region pixels that do not belong to the green pixel sets
			sumi += (*it)[0];
			sumj += (*it)[1];
		}
		sumi /= n;
		sumj /= n;

		this->centroid[0] = sumi;
		this->centroid[1] = sumj;
	}


}

/**
 * @return class of one of the 5 phagocythosis types 
 */
cell_class::pathogen RegionP::getClassP(){
	return this->klassP;
}

// TODO das so machen : double area = fabs(contourArea(cv::Mat(contours[i])));
// if( fabs( area ) >= (double) p){
//     cv::drawContours(dst, contours, i, cv::Scalar(255), -1);
// }
int RegionP::getAreaGray(){ 
	return this->region_pixels.size() - this->green_region_pixels.size(); 
}

/**
 * @return region centroid as cv::Scalar
 */
cv::Scalar RegionP::getCentroidScalar(){
	return cv::Scalar(centroid[0], centroid[1]);
}

/**
 * @return region centroid as cv::Point
 */
cv::Point RegionP::getCentroidPoint(){
	return cv::Point(centroid[1], centroid[0]);
}

void RegionP::getDeadFungiIds(std::vector<int> & ids){
	ids = this->deadFungiIDs;
}


int RegionP::getIState(){
	return this->I_state;
}

bool RegionP::getFCState(){
	return this->FC;
}

double RegionP::getFungalRatioGreen(){
	return (double) this->green_region_pixels.size() / (double) this->getArea();
}

void RegionP::getFungiIds(std::vector<int> & ids){
	ids = this->fungiIDs;
}

int RegionP::getNumberOfPhagocytoses(){
	return this->phagocytosisIDs.size();
}

void RegionP::getPhagocytosisIds(std::vector<int> & ids){
	ids = this->phagocytosisIDs;
}

bool RegionP::hasFungiID(int id){
	if(std::find(this->fungiIDs.begin(), this->fungiIDs.end(), id) == this->fungiIDs.end()){
		return false;
	}
	return true;
}

bool RegionP::hasPhagocytozedFungus(int id){
	if(std::find(this->phagocytosisIDs.begin(), this->phagocytosisIDs.end(), id) == this->phagocytosisIDs.end()){
		return false;
	}
	return true; //else
}

void RegionP::removeFungiId(int id){
	std::vector<int>::iterator it = std::find(this->fungiIDs.begin(), this->fungiIDs.end(), id);
	if(it != this->fungiIDs.end()){
		this->fungiIDs.erase(it);
	}
}

void RegionP::removeDeadFungiIds(){
	this->deadFungiIDs.clear();
}

void RegionP::removeFungiIds(){
	this->fungiIDs.clear();
}


void RegionP::removeLastPhID(){
	if(this->phagocytosisIDs.size() > 0){
		this->phagocytosisIDs.pop_back();
	}
}

void RegionP::removePhIDs(){
	if(this->phagocytosisIDs.size() > 0){
//		cout<< this->phagocytosisIDs.size();
		this->phagocytosisIDs.resize(0);
	}
}

void RegionP::setClassP(cell_class::pathogen c){
	this->klassP = c;
}

void RegionP::setIState(int state){
//	cout << "set istate " << state << endl;

	this->I_state=state;
}

/**
 * This function sets the class of a region according to: 0 -> noise, 1 -> single cell, 2 -> cell cluster
 *
 * @param k class
 */
void RegionP::setKlass(cell_class::immune k){
	this->klass = k;
}

void RegionP::setFCState(bool fc){
	this->FC=fc;
}