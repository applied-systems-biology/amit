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

#include "Calcs.h"

namespace Calcs 
{

    /**
     * gets min and max values of the two input vectors
     * time complexity: O(#region pixels)
     */
    void minMaxVec(std::vector<cv::Vec2i> &v, std::vector<cv::Vec2i> &w, cv::Vec2i &min, cv::Vec2i &max) {

        std::vector<cv::Vec2i>::iterator it_v, it_w;

        // calculate minimum and maximum x- and y-coordinates
        int minx = INT_MAX;
        int miny = INT_MAX;
        int maxx = 0, maxy = 0;

        // check first vector
        for(it_v = v.begin(); it_v != v.end(); it_v++){
            int x = it_v->val[0];
            int y = it_v->val[1];

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

        // check second vector
        for(it_w = w.begin(); it_w != w.end(); it_w++){
            int x = it_w->val[0];
            int y = it_w->val[1];

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

        min[0] = minx;
        min[1] = miny;
        max[0] = maxx;
        max[1] = maxy;

    }

    /**
     * This function computes the mean of a set of double values.
     *
     * @param v set of double values
     *
     * @return mean
     */
    double computeMean(const std::vector<double> &v){
        std::vector<double>::const_iterator it;
        int n = v.size();
        double sum = 0.0;

        for(it = v.begin(); it != v.end(); it++){
            sum += (*it);
        }

        sum /= n;

        return sum;
    }

    /**
     * This function computes the variance of a set of double values.
     *
     * @param v set of double values
     *
     * @return variance
     */
    double computeVariance(const std::vector<double> &v){
        double var = 0.0;
        double m = Calcs::computeMean(v);

        for(std::vector<double>::const_iterator it = v.begin(); it != v.end(); it++){
            var += pow(((*it) - m), 2.0);
        }

        var /= (v.size()-1);
    
        return var;
    }

    /**
     * This function computes the standard deviation of a set of double values.
     *
     * @param v set of double values
     *
     * @return standard deviation
     */
    double computeSd(const std::vector<double> &v){
        return sqrt(Calcs::computeVariance(v));
    }

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

        if(it == list.end()){
            return false;
        }
        else{
            return true;
        }
        
    }

} // Calcs