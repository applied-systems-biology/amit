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

#include "Region.h"
#include "Calcs.h"


namespace cell_class 
{
    
    std::string getCellClass(immune c){
        switch(c)
        {
            case immune::NOISE      : return "N";
            case immune::SINGLE     : return "S";
            case immune::CLUSTER    : return "C";
            default                 : return "ERROR";
        }
    };

    immune setCellType(int c){
        switch(c)
        {
            case 0  : return immune::NOISE;
            case 1  : return immune::SINGLE;
            case 2  : return immune::CLUSTER;
            default : return immune::SINGLE;
        }
    }

	pathogen setCellTypeP(int c){
		switch(c)
        {
            case 1  : return pathogen::NON_INTERACTING;
            case 2  : return pathogen::IMMUNE_TOUCHING_PATHOGEN;
            case 3  : return pathogen::PHAGOCYTOSED;
			case 4  : return pathogen::PHAGOCYTOSED_INTERACTION_IMMUNE;
            case 5  : return pathogen::PHAGOCYTOSED_INTERACTION_PATHOGEN;
			case 6  : return pathogen::UNDEFINED;
            default : return pathogen::NON_INTERACTING;
        }
	}

};


/**
 * This function checks if a given pixel is present in a vector of pixels.
 *
 * @param v pixels
 * @param list set of pixels
 *
 * @return true, if the pixels is found in the set of pixels
 */
bool exists(cv::Vec2i &v, std::vector<cv::Vec2i> &list){
	std::vector<cv::Vec2i>::iterator it;
	it = find(list.begin(), list.end(), v);

	if(it == list.end())
		return false;
    else
        return true;
}

///// class member functions /////

/**
 * This function returns the class of the region of an immune cell.
 *
 * @return class
 */
cell_class::immune Region::getClass(){
	return this->klass;
}

/**
 * @return ID
 */
int Region::getId(){
	return this->ID;
}

/**
 * @return region area // TODO das so machen : double area = fabs(contourArea(cv::Mat(contours[i])));

            if( fabs( area ) >= (double) p){
                cv::drawContours(dst, contours, i, cv::Scalar(255), -1);
            }
 */
int Region::getArea(){
	return (int) this->region_pixels.size();
}

/**
 * This function returns the minimal enclosing image.
 *
 * @return minimal enclosing image
 */
cv::Mat Region::getImage(){
	this->createImage();
	return this->im;
}

/**
 * @return region centroid as cv::Scalar
 */
cv::Scalar Region::getCentroidScalar(){
	return cv::Scalar(centroid[0], centroid[1]);
}

/**
 * @return region centroid as cv::Point
 */
cv::Point Region::getCentroidPoint(){
	return cv::Point(centroid[1], centroid[0]);
}

/**
 * This function returns minimal and maximal x- and y-coordinates of the region.
 *
 * @param minx minimal x-coordinate
 * @param miny minimal y-coordinate
 * @param maxx maximal x-coordinate
 * @param maxy maximal y-coordinate
 */
void Region::getMinMax(int &minx, int &miny, int &maxx, int &maxy){
	minx = INT_MAX;
	miny = INT_MAX;
	maxx = 0;
	maxy = 0;

	std::vector<cv::Vec2i>::iterator vecit;
	for(vecit = this->contour_pixels.begin(); vecit != this->contour_pixels.end(); vecit++){
		int x = vecit->val[0];
		int y = vecit->val[1];

		if(x < minx)
			minx = x;
		if(x > maxx)
			maxx = x;
		if(y < miny)
			miny = y;
		if(y > maxy)
			maxy = y;
	}
}

/**
 * This function computes the minimum enclosing ellipse and returns its principal axes.
 *
 * @param s principal axes of minimum enclosing ellipse
 */
void Region::getPrincipalAxesMinEnclosingEllipse(cv::Size2f &s){

	// contour size of at least 5 pixels is necessary to compute enclosing ellipse
	if(this->contour_pixels.size() < 5){
		s.width = 0;
		s.height = 0;
	}
	else{

		std::vector<cv::Point> points;
		for(std::vector<cv::Vec2i>::iterator it = this->contour_pixels.begin(); it != this->contour_pixels.end(); it++){
			cv::Point p(it->val[1], it->val[0]);
			points.push_back(p);
		}

		cv::RotatedRect rect = cv::fitEllipse(points);
		s = rect.size;
	}
}

/**
 * return true, if region is interacting (in contact) with another region
 */
bool Region::isInteracting(){
	if(this->interactionIDs.size() > 0){
		return true;
	}
	return false;
}

/**
 * This function gives the IDs of interacting cells.
 *
 * This is a call-by-reference, changes made to the parameter affect the passed argument!
 *
 * @param ids set of IDs is returned
 */
void Region::getInteractionIds(std::vector<int> &ids){
	ids = this->interactionIDs;
}

/**
 * @return minor axe of ellipsoid that fits the region
 */
float Region::getMin(){
	return this->minAxe;
}

/**
 * @return major axe of ellipsoid that fits the region
 */
float Region::getMaj(){
	return this->majAxe;
}

/**
 * @return the status the region's flattness
 */
bool Region::getFlat(){
	return this->flat;
}

/**
 * This function adds a vector of pixels to the set of region pixels.
 *
 * @param v set of pixels
 */
void Region::addRegionPixels(const std::vector<cv::Vec2i> &v){
	std::vector<cv::Vec2i>::const_iterator it;
	for(it = v.begin(); it != v.end(); it++){
		this->region_pixels.push_back((*it));
	}
}

/**
 * This function adds a vector of pixels to the set of contour pixels.
 *
 * @param v set of pixels
 */
void Region::addContourPixels(const std::vector<cv::Vec2i> &v){
	std::vector<cv::Vec2i>::const_iterator it;
	for(it = v.begin(); it != v.end(); it++){
		this->contour_pixels.push_back((*it));
	}
}

/**
 * This function adds a single pixel to the set of region pixels.
 *
 * @param p pixel in cv::Vec2i format
 */
void Region::addRegionPixel(const cv::Vec2i &p){
	this->region_pixels.push_back(p);
}

/**
 * This function includes the ID of an interacting cell.
 * IDs cannot be added more than once.
 * The list of IDs is sorted in numerically increasing mode.
 *
 * @param id ID
 */
void Region::addInteractionId(const int &id){
	if(std::find(this->interactionIDs.begin(), this->interactionIDs.end(), id) == this->interactionIDs.end()){
		this->interactionIDs.push_back(id);
		sort(this->interactionIDs.begin(), this->interactionIDs.end());
	}
}

/**
 * This function checks if a pixel has a neighbor in the Moore-neighborhood within a set of pixels.
 *
 * @param v pixel
 * @param list set of pixels
 *
 * return true, if pixel has a neighbor in the set of pixels
 */
bool Region::hasNeighbor(cv::Vec2i &v, std::vector<cv::Vec2i> &list){
	//(x-1|y)
	cv::Vec2i v1(v[0]-1, v[1]);
	if(!exists(v1, list)){
		return true;
	}

	//(x-1|y-1)
	cv::Vec2i v2(v[0]-1, v[1]-1);
	if(!exists(v2, list)){
		return true;
	}

	//(x|y-1)
	cv::Vec2i v3(v[0], v[1]-1);
	if(!exists(v3, list)){
		return true;
	}

	//(x+1|y-1)
	cv::Vec2i v4(v[0]+1, v[1]-1);
	if(!exists(v4, list)){
		return true;
	}

	//(x+1|y)
	cv::Vec2i v5(v[0]+1, v[1]);
	if(!exists(v5, list)){
		return true;
	}

	//(x+1|y+1)
	cv::Vec2i v6(v[0]+1, v[1]+1);
	if(!exists(v6,list)){
		return true;
	}

	//(x|y+1)
	cv::Vec2i v7(v[0], v[1]+1);
	if(!exists(v7, list)){
		return true;
	}


	//(x-1|y+1)
	cv::Vec2i v8(v[0]-1, v[1]+1);
	if(!exists(v8, list)){
		return true;
	}


	return false;
}

/**
 * This function computes the number of pixels that are shared by two regions.
 *
 * @param r first region
 * @param p second region
 *
 * @return number of shared pixels
 */
int overlapN(Region &r, Region &p){
	std::vector<cv::Vec2i>::iterator it_r, it_p;

	cv::Scalar c1 = r.getCentroidScalar();
	cv::Scalar c2 = p.getCentroidScalar();

	int x_r = c1[0];
	int x_p = c2[0];
	int y_r = c1[1];
	int y_p = c2[1];

	if(sqrt(pow(x_r - x_p, 2.0) + pow(y_r - y_p, 2.0)) > 500){
		return 0;
	}
	else {
		cv::Vec2i min, max;
		Calcs::minMaxVec(r.region_pixels, p.region_pixels, min, max);

		// create images that include both regions but not more
		int imrows = max[0] - min[0] + 1;
		int imcols = max[1] - min[1] + 1;

		cv::Mat tempr(imrows, imcols, CV_8UC1, cv::Scalar(0));
		cv::Mat tempp(imrows, imcols, CV_8UC1, cv::Scalar(0));

		// create images with only one region (r and p)
		for(it_r = r.region_pixels.begin(); it_r != r.region_pixels.end(); it_r++){
			int x = it_r->val[0] - min[0];
			int y = it_r->val[1] - min[1];
			tempr.at<uchar>(x,y) = 255;
		}

		for(it_p = p.region_pixels.begin(); it_p != p.region_pixels.end(); it_p++){
			int x = it_p->val[0] - min[0];
			int y = it_p->val[1] - min[1];
			tempp.at<uchar>(x,y) = 255;
		}

		cv::Mat overlap(imrows, imcols, CV_8UC1);
	
		cv::bitwise_and(tempr, tempp, overlap);

		int s = cv::sum(overlap)[0]/255;

		return s;
	}
}

/**
 * This function checks if two regions share a common set of pixels.
 *
 * @param r first region
 * @param p second region
 *
 * @return true, if an overlap exists, else, false
 */
//bool overlap(Region &r, Region &p){
int overlap(Region &r, Region &p){
	return overlap(&r, &p);
}

/**
 * This function checks if two regions share a common set of pixels.
 *
 * @param r first region
 * @param p second region
 *
 * @return true, if an overlap exists, else, false
 */
//bool overlap(Region *r, Region *p){
int overlap(Region *r, Region *p){

	std::vector<cv::Vec2i>::iterator it_r, it_p;

	cv::Scalar c1 = r->getCentroidScalar();
	cv::Scalar c2 = p->getCentroidScalar();

	int x_r = c1[0];
	int x_p = c2[0];
	int y_r = c1[1];
	int y_p = c2[1];

	// CAUTION / WARNING: hard threshold by euclidian distance of 500 px
	if(sqrt(pow(x_r - x_p, 2.0) + pow(y_r - y_p, 2.0)) > 500){
		return false;
	}

	cv::Vec2i min, max;
	Calcs::minMaxVec(r->region_pixels, p->region_pixels, min, max);

	// create images that include both regions but not more
	int imrows = max[0] - min[0] + 1;
	int imcols = max[1] - min[1] + 1;

	cv::Mat tempr(imrows, imcols, CV_8UC1, cv::Scalar(0));
	cv::Mat tempp(imrows, imcols, CV_8UC1, cv::Scalar(0));

	// create images with only one region (r and p)
	for(it_r = r->region_pixels.begin(); it_r != r->region_pixels.end(); it_r++){
		int x = it_r->val[0] - min[0];
		int y = it_r->val[1] - min[1];
		tempr.at<uchar>(x,y) = 255;
	}

	for(it_p = p->region_pixels.begin(); it_p != p->region_pixels.end(); it_p++){
		int x = it_p->val[0] - min[0];
		int y = it_p->val[1] - min[1];
		tempp.at<uchar>(x,y) = 255;
	}

	cv::Mat overlap(imrows, imcols, CV_8UC1);

	cv::bitwise_and(tempr, tempp, overlap);
	
	int s = cv::sum(overlap)[0];

	/// regions overlap with the sum over all pixels intensities
	if (s > 0) {
//        return true;
        return s;
    }
	else {
        return 0; //false;
	}

}

/**
 * This function sets the class of a region according to: 0 -> noise, 1 -> single cell, 2 -> cell cluster
 *
 * @param k class
 */
void Region::setClass(cell_class::immune c){
	this->klass = c;
}

/**
 * This function sets the ID of the region.
 *
 * @param ID ID
 */
void Region::setId(const int &id){
	this->ID = id;
}

/**
 * This function sets the centroid of the region
 *
 * @param x x-coordinate
 * @param y y-coordinate
 */
void Region::setCentroid(const double &x, const double &y){
	this->centroid.val[0] = x;
	this->centroid.val[1] = y;
}

/**
 * This function sets the flat status of the region
 *
 * @param f flat flag
 */
void Region::setFlat(bool f){
	this->flat = f;
}

/**
 * This function sets the major and minor axes of ellipsoid that fits the region
 *
 * @param min minor axe
 * @param maj major axe
 */
void Region::setMinMajAxe(const float min, const float maj){
	this->minAxe = min;
	this->majAxe = maj;
}

/**
 * This function computes the contour (Moore-neighborhood) of the region.
 */
void Region::computeContour(){
	std::vector<cv::Vec2i>::iterator it;
	this->contour_pixels.resize(0);

	for(it = this->region_pixels.begin(); it != this->region_pixels.end(); it++){
		if(Region::hasNeighbor(*it, this->region_pixels)){
			this->contour_pixels.push_back(*it);
		}
	}
}

/**
 * This function computes the center of mass from all pixels stored in the vector of region pixels.
 */
void Region::computeCentroid(){
	std::vector<cv::Vec2i>::iterator it;
	int n =  this->getArea();

	double sumi = 0.0;
    double sumj = 0.0;

	for(it = this->region_pixels.begin(); it != this->region_pixels.end(); it++){
		sumi += it->val[0];
		sumj += it->val[1];
	}
	sumi /= n;
	sumj /= n;

	this->centroid[0] = sumi;
	this->centroid[1] = sumj;
}

/**
 * returns the distance between the centroids of two regions
 */
double computeDistance(Region &p, Region &q){
	cv::Scalar c1 = p.getCentroidScalar();
	cv::Scalar c2 = q.getCentroidScalar();

	double dist = sqrt(pow((double)(c1[0] - c2[0]),2.) + pow((double)(c1[1] - c2[1]),2.));
	return dist;
}

/**
 * returns the distance between the centroids of two regions
 */
double computeDistance(Region *p, Region *q){
	cv::Scalar c1 = p->getCentroidScalar();
	cv::Scalar c2 = q->getCentroidScalar();

	double dist = sqrt(pow((double)(c1[0] - c2[0]),2.) + pow((double)(c1[1] - c2[1]),2.));
	return dist;
}

/**
 * This function computes the minimal enclosing rectangle for the region and saves this as an image including the region as foreground.
 * TODO: do this with cv::boundingRect
 */
void Region::createImage(){
	std::vector<cv::Vec2i>::iterator vecit;

	int minx = INT_MAX;
	int miny
     = INT_MAX;
	int maxx = 0, maxy = 0;

	for(vecit = this->contour_pixels.begin(); vecit != this->contour_pixels.end(); vecit++){
		int x = vecit->val[0];
		int y = vecit->val[1];

		if(x < minx){
			minx = x;
		}
		if(x > maxx){
			maxx = x;
		}
		if(y < miny){
			miny = y;
		}
		if(y > maxy){
			maxy = y;
		}
	}

	cv::Mat tmp(maxx-minx+1, maxy-miny+1, CV_8UC1, cv::Scalar(0));
	this->im = tmp;

	for(vecit =  this->region_pixels.begin(); vecit != this->region_pixels.end(); vecit++){
		cv::Point p(vecit->val[1]-miny, vecit->val[0]-minx);
		im.at<uchar>(p) = 255;
	}

}

/**
 * This function changes the ID id of interacting cells to a new ID id_new
 *
 * @param id old ID
 * @param id_new new ID
 */
void Region::changeInteractionID(int id, int id_new){
	std::replace(this->interactionIDs.begin(), this->interactionIDs.end(), id, id_new);
}


