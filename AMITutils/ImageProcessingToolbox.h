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

/*
This namespace contains several functions of a script language image processing toolbox that were transferred to C++ using the OpenCV library
*/

#pragma once

#include <vector>
#include <string>
#include <optional>
#include <opencv2/opencv.hpp>


namespace ImageProcessingToolbox
{
    struct tabulate
    {
        std::vector<int> value;
        std::vector<int> count;
        std::vector<float> percent;

        // create absolute / relative frequency map
        tabulate(const std::vector<int> &values) {
            std::map<int,int> frequency_count;
            for (size_t i = 0; i < values.size(); i++)
                frequency_count[values[i]]++;            

            // compute all vectors
            for(std::map<int,int>:: iterator it = frequency_count.begin(); it != frequency_count.end(); it++)
            {
                this->value.push_back( it->first );
                this->count.push_back( it->second );
                this->percent.push_back( (float) it->second / values.size() );            
            }
        }

    };

    // functions

    std::string type2str(int type);
    cv::Mat createKernel(const std::string kernelType);
    std::vector<float> unique(const cv::Mat& input, bool sort = false);

    void imgradient(const cv::Mat &src, cv::Mat &magnitude);
    void imoverlay(const cv::Mat &src, const cv::Mat &mask, cv::Mat &dst, const std::string color);
    void imfill(const cv::Mat &src, cv::Mat &dst);
    void imclearborder(const cv::Mat &src, cv::Mat &dst, const int &conn = 4);
    void imerode(const cv::Mat &src, cv::Mat &dst, const cv::Mat &kernel, const int &iterations = 1);
    void imdilate(const cv::Mat &src, cv::Mat &dst, const cv::Mat &kernel, const int &iterations = 1);
    void imclose(const cv::Mat &src, cv::Mat &dst, const cv::Mat &kernel, const int &iterations = 1);
    void imopen(const cv::Mat &src, cv::Mat &dst, const cv::Mat &kernel, const int &iterations = 1);
    void imadjust(const cv::Mat1b &src, cv::Mat1b &dst, int tol = 1, cv::Vec2i in = cv::Vec2i(0, 255), cv::Vec2i out = cv::Vec2i(0, 255));
    void imbinarize(const cv::Mat &src, cv::Mat1b &dst, const double &globThreshold=-1.0, const bool &foregroundBright=true, const int &offset=255, const float &sensitivity=13.);
    void imsubtract(const cv::Mat &minuating, const cv::Mat &subtracting, cv::Mat &dst);
    void imadd(const cv::Mat &summand1, const cv::Mat &summand2, cv::Mat &dst);
    void bwlabel(const cv::Mat1b &src, cv::Mat &dst, int &nLabels, int connectivity);
    void bwmorph(const cv::Mat &src, cv::Mat &dst, const std::string operation);
    void bwareaopen(const cv::Mat &src, cv::Mat &dst, const int &p);
    void regionprops(const cv::Mat &src, std::vector<std::vector<cv::Point> > &object_coordinates, const std::string &method);
    void adapthisteq(const cv::Mat &src, cv::Mat &dst, const double clipLimit=4.0, const std::string distribution="uniform");
    void stdfilt(const cv::Mat &src, cv::Mat &dst);
    
    // needed licence
    void bwlookup( const cv::Mat & in, cv::Mat & out, const cv::Mat & lut, int bordertype=cv::BORDER_CONSTANT, cv::Scalar px = cv::Scalar(0) ); 
    void wiener2(const cv::Mat &src, cv::Mat &dst, const cv::Size &kSize, double &noiseVariance);

}