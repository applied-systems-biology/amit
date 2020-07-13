/*  
Copyright by Stefanie Dietrich 

Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo
Figge
https://www.leibniz-hki.de/en/applied-systems-biology.html
HKI-Center for Systems Biology of Infection
Leibniz Institute for Natural Product Research and Infection Biology -
Hans Knöll Insitute (HKI)
Adolf-Reichwein-Straße 23, 07745 Jena, Germany 

This code is licensed under BSD 2-Clause
See the LICENSE file provided with this code for the full license.
 */

#pragma once

#include <string>
#include <vector>
#include <opencv2/opencv.hpp>


namespace singleCells {

    cv::Mat processSingle(const cv::Mat &src, const bool &FLAG_DEBUG);

    void findContourPoints(cv::Mat bg_img, std::vector<std::vector<cv::Point>> contours, std::vector<cv::Vec4i> hierarchy, int min_area, int min_dst, int n1, int m, std::vector<std::vector<cv::Point>> &cont_array, std::vector<std::vector<cv::Point>> &midpoint_array);
    void getLinePoints(int j, cv::Point pt1, cv::Point pt2, std::vector<cv::Point> &linePoints);
    int myRounding(double val);

    void findCutLines(std::vector<std::vector<cv::Point>> cont_array, std::vector<std::vector<cv::Point>> midpoint_array, std::vector<std::vector<cv::Point>> contours, std::vector<cv::Vec4i> hierarchy, int min_width, int max_width, int min_length, int max_length, cv::Mat bg_img, std::vector<std::vector<std::vector<cv::Point>>> &cut_lines_array);
    void getRect(cv::Point cp, cv::Point mp, double width, double length, std::vector<cv::Point> &rect);
    cv::Point findPoint(int ix, cv::Point cp, std::vector<cv::Point> rect, std::vector<cv::Point> set1, bool& found );

    void postprocessCutLines(std::vector<std::vector<std::vector<cv::Point>>> &cut_lines_array);
    std::vector<cv::Point> getPoints(cv::Point point, std::vector<std::vector<cv::Point>> lines, bool& triangle, std::vector<int>& num);
    void hashSearchPoint(cv::Point  hashtable[], cv::Point point, long &pos, bool &notcontained, int &nn);
    long hashvalue(cv::Point punkt);
    std::vector<int> decimal2binary(long int z);

    void singleCells_segmentation(const std::vector<cv::Mat> &images, const int &N_THREADS, const std::string &OUT_BINARY, const bool &FLAG_DEBUG);

} // singleCells