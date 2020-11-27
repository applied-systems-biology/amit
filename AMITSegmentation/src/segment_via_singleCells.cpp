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

#include "segment_via_singleCells.h"
#include <iostream>
// #include <opencv2/highgui.hpp>
// #include <opencv2/imgproc.hpp>
// #include <opencv2/imgcodecs.hpp>
// #include "opencv2/core.hpp"
#include "IOputs.h"
#include "ImageProcessingToolbox.h"

namespace IPT = ImageProcessingToolbox;


namespace singleCells
{

    cv::Mat processSingle(const cv::Mat &src, const bool &FLAG_DEBUG){

        /***** preprocessing *****/
        
        /// *** Minimum - filtering *** ///
        cv::Mat img_min, img_min_cvt;

        cv::Mat min_kernel = cv::getStructuringElement( cv::MORPH_RECT, cv::Size(3,3) ); 
        cv::erode(src, img_min, min_kernel);

        /// convert to 1-channel image (grayscaled)
        cv::cvtColor(img_min, img_min_cvt, cv::COLOR_BGR2GRAY);

        // TODO: add gamma calculation

        /// *** Median - blur filtering *** ///
        cv::Mat img_med;

        cv::medianBlur(img_min_cvt, img_med, 3);

        // TODO: add gamma calculation

        /// *** Opening *** ///
        cv::Mat img_open;

        cv::Mat open_kernel = cv::getStructuringElement( cv::MORPH_ELLIPSE, cv::Size(5,5) ); 
        cv::morphologyEx(img_med, img_open, cv::MORPH_OPEN, open_kernel);

        /// *** Normalize *** ///
        cv::Mat img_norm;

        double  minVal,  maxVal;
        cv::minMaxLoc(img_min_cvt, &minVal, &maxVal); 
        img_open.convertTo(img_norm, CV_8UC1, (1.0 / maxVal)*255. , 0. );
        
        /// *** Otsu - thresholding with preceding Gaussian blur *** ///
        cv::Mat img_blur, img_blur_norm, img_blur_otsu;

        cv::GaussianBlur(img_norm, img_blur, cv::Size(3,3), 2);
        /// normalize again
        cv::minMaxLoc(img_blur, &minVal, &maxVal); 
        img_blur.convertTo(img_blur_norm, CV_8UC1, (1.0 / maxVal)*255. , 0. );

        cv::threshold(img_blur_norm, img_blur_otsu, cv::THRESH_OTSU, 255, cv::THRESH_BINARY);

        /// ### Convexity and separation nach Farhan 2013 ### ///
            
        const int n = 20;            // length of tendon for the search of convexity points
        const int m = 6;             // TODO übersetzen     # abstand der endpunkte der sehne zum convexity punkt
        const int min_dst = 1;       // TODO übersetzen     # minimaler abstand des convexity punkt zu seiner sehne
        const int min_area = 60;     // TODO übersetzen     # mindest clump groesse, um bei der segmentierung betrachtet zu werden
        
        const int min_width = 6;     // TODO übersetzen     # minimale Breite des Suchrechtecks
        const int min_length = 30;   // TODO übersetzen     # minimale Laenge des Suchrechtecks
        const int max_width = 30;    // TODO übersetzen     # maximale Breite des Suchrechtecks
        const int max_length = 40;   // TODO übersetzen      # maximale Laenge des Suchrechtecks


        std::vector<std::vector<cv::Point>> contours;
        std::vector<cv::Vec4i> hierarchy;

        cv::findContours(img_blur_otsu, contours, hierarchy, cv::RETR_CCOMP, cv::CHAIN_APPROX_NONE);
       
        
        /***** Find Contour Points *****/
        std::vector<std::vector<cv::Point>> cont_array;
	    std::vector<std::vector<cv::Point>> midpoint_array;

        cv::Mat bg_img;
        img_blur_otsu.copyTo(bg_img);

        singleCells::findContourPoints(bg_img, contours, hierarchy, min_area, min_dst, n, m, cont_array, midpoint_array);
        
        std::cout << "cont_array-size: " << cont_array.size() << std::endl;
        std::cout << "midpoint_array-size: " << midpoint_array.size() << std::endl;

        /******** Find Cut-lines *******/
        std::vector<std::vector<std::vector<cv::Point>>> cut_lines_array;

        singleCells::findCutLines(cont_array, midpoint_array, contours, hierarchy, min_width, max_width, min_length, max_length, bg_img, cut_lines_array);

        std::cout << "cut_lines_array-size: " << cut_lines_array.size() << std::endl;

        /**** Postprocess cut lines ****/
        postprocessCutLines(cut_lines_array);
	    
        std::cout << "post processed cut lines: " << cut_lines_array.size() << std::endl;
        
        /********** Draw Lines *********/
        cv::Mat img_singlecells = bg_img.clone();

        for (int i = 0; i < (int)cut_lines_array.size(); i++){
            for (int j = 0; j < (int)cut_lines_array[i].size(); j++){
                cv::line(img_singlecells, cut_lines_array[i][j][0], cut_lines_array[i][j][1], cv::Scalar(0), 2, 4, 0);
            }
        }
        
        std::cout << "Draw cut lines" << std::endl;

        /******* Postprocessing 2 ******/

        // TODO: hier evtl convertTo function nutzen um zu normalisieren
        cv::normalize(img_singlecells, img_singlecells, 0, 255, cv::NORM_MINMAX, CV_8U);
        
        contours.clear();
        hierarchy.clear();
        cv::findContours(img_singlecells.clone(), contours, hierarchy, cv::RETR_CCOMP, cv::CHAIN_APPROX_SIMPLE);
        
        int num_cells = (int)contours.size();
        double area;
        
        std::vector<cv::Point> cnt;
        for (int i = 0; i < (int)contours.size(); i++){
            cnt = contours[i];
            area = cv::contourArea(cnt);
            if (area < 30){
                cv::drawContours(img_singlecells, contours, i, cv::Scalar(0), cv::FILLED, 8);
                num_cells -= 1;
            }
        }

        std::cout << "Second post processing done" << std::endl;





        /// apply gamma (0.3) on intensities
        // cv::Scalar gamma = cv::Scalar(0.3); 
        // img_min_gamma = cv::exp()

        // double a[] = { 1,2,3,4,5,6,7,8,9 };
        // cv::Mat m0(3, 3, CV_64FC1, a);
        // cv::Mat m1(3, 3, CV_64FC1);
        // cv::Mat m2(3, 3, CV_64FC1);
        // std::cout << "Matrix m0 : \n" << m0 << "\n";
        // cv::exp(m0, m1);
        // std::cout << "Matrix m1 : \n" << m1 << "\n";

        std::cout << "img_blur_otsu-contours: " << contours.size() << std::endl;

        
        cv::minMaxLoc(img_blur, &minVal, &maxVal); 
        std::cout << "\n img_blur  : min - max: " << minVal << " - " << maxVal << std::endl;
        
        cv::minMaxLoc(img_blur_otsu, &minVal, &maxVal); 
        std::cout << "\n img_blur_otsu  : min - max: " << minVal << " - " << maxVal << std::endl;
        
        // std::cout << "img_min-type: " << IPT::type2str(img_min.type()) << std::endl;
        // std::cout << "img_min_cvt-type: " << IPT::type2str(img_min_cvt.type()) << std::endl;
        // std::cout << "img_med-type: " << IPT::type2str(img_med.type()) << std::endl;
        // std::cout << "img_open-type: " << IPT::type2str(img_open.type()) << std::endl;
        // std::cout << "img_norm-type: " << IPT::type2str(img_norm.type()) << std::endl;
        // std::cout << "img_blur-type: " << IPT::type2str(img_blur.type()) << std::endl;
        std::cout << "img_blur_norm-type: " << IPT::type2str(img_blur_norm.type()) << std::endl;
        std::cout << "img_blur_otsu-type: " << IPT::type2str(img_blur_otsu.type()) << std::endl;


        // cv::imshow("src", src);
        // cv::imshow("img_min", img_min);
        // cv::imshow("img_min_cvt", img_min_cvt);
        // cv::imshow("img_med", img_med);
        // cv::imshow("img_open", img_open);
        // cv::imshow("img_norm", img_norm);
        // cv::imshow("img_blur", img_blur);
        cv::imshow("img_blur_norm", img_blur_norm);
        cv::imshow("img_blur_otsu", img_blur_otsu);
        cv::imshow("img_singlecells", img_singlecells);
		cv::waitKey(0);

        return img_singlecells;

    }

    /**
     * finds the concavity points in clumps
     * bg_img = b/w image of clumps
     * contours = contours of the clumps
     * hierarchy = hierarchy of the contours
     * min_area = minimal clumps size
     * min_dst = minimal distance from concavity point to chord (minimal depth of convexity defect)
     * n1 = length of the contour segment for concavity point search (smaller if inside whole)
     * m = half of the chord over a concaviyt point
    */
    void findContourPoints(cv::Mat bg_img, std::vector<std::vector<cv::Point>> contours, std::vector<cv::Vec4i> hierarchy, int min_area, int min_dst, int n1, int m, std::vector<std::vector<cv::Point>> &cont_array, std::vector<std::vector<cv::Point>> &midpoint_array){
        
        cv::normalize(bg_img, bg_img, 0, 255, cv::NORM_MINMAX);// TODO: evtl mit convertTo function machen

        for (int i=0; i < (int)contours.size(); i++){

            std::cout << "contour number: " << i << std::endl;

            std::vector<cv::Point> cnt = contours[i];

            double area = cv::contourArea(cnt);
            int n;
            
            std::vector<cv::Point> con_points;
            std::vector<cv::Point> midpoints;

            if (hierarchy[i][3] != -1){		// if inside a hole
                n = 10;
            }
            else {
                n = n1;
            } 

            if (area > min_area){

                std::cout << "Area is bigger than min_area." << std::endl;

                // cnt2 will have beginning-values mirrored at the end
                std::vector <cv::Point> cnt2;		
                
                for (int k = 0; k < (int)cnt.size(); k++){ // TODO: mit size_t ersetzen
                    cnt2.push_back(cnt[k]);
                }
                for (int k = (int)cnt.size(); k < (int)cnt.size()+m+n; k++){
                    cnt2.push_back(cnt[k-(int)cnt.size()]);
                }

                int j = 0;

                std::cout << "created longer contour vector. \n"
                        "go through this contour vector next:" << std::endl;

                while (j < (int)cnt.size()-(n/2)){
                    cv::Point pt1 = cnt2[j];
                    cv::Point pt2 = cnt2[j+n];

                    // contour points between the two points
                    std::vector<cv::Point> curve(n+1);	
                    
                    for (int k = j; k < j+n+1; k++){
                        curve[k-j] = cnt2[k];
                    }

                    std::cout << "- get line between next two contour points..." << std::endl;
                    std::vector<cv::Point> line;
                    singleCells::getLinePoints(j,pt1, pt2, line);

                    std::cout << "- check if line crosses background..." << std::endl;
                    int summe = 0;
                    for (int k = 0; k < (int)line.size(); k++){
                        summe += bg_img.at<uchar>(line[k]);
                    }
                    
                    std::cout << "- if line crosses background, check all points in the background..." << std::endl;
                    if ((double)summe/255. <= (int)line.size()-2){
                        
                        std::cout << "- - - gather all points in background." << std::endl;
                        std::vector<int> points;
                        
                        // get all points in BG
                        for (int k = 0; k < (int)line.size(); k++){	
                            if (bg_img.at<uchar>(line[k]) == 0){
                                points.push_back(k);
                            }
                        }

                        cv::Point p1 = line[points[0]-1];
                        cv::Point p2 = line[points[points.size()-1]+1];

                        double sum1 = std::numeric_limits<double>::infinity();
                        double sum2 = std::numeric_limits<double>::infinity();
                        int c1 = 0;
                        int c2 = 0;

                        std::cout << "- - - 5" << std::endl;

                        for (int k = 0; k < (int)curve.size(); k++){
                            double temp = abs(curve[k].x-p1.x)+ abs(curve[k].y-p1.y);
                            if (temp < sum1){
                                c1 = k;
                                sum1 = temp;
                            }
                            temp = abs(curve[k].x-p2.x) + abs(curve[k].y-p2.y);
                            if (temp  < sum2){
                                c2 = k;
                                sum2 = temp;
                            }
                        }

                        std::cout << "- - - 6" << std::endl;
                        c2 = curve.size()-1-c2;

                        std::vector<cv::Point> curve2(curve.begin()+c1, curve.end()-c2);
                        
                        double win_dst = 0;
                        cv::Point winner;
                        cv::Point winner_middle;
                        int pos = 0;
                        std::cout << "- - - go through all points in curve2 and find the winner" << std::endl;
                        
                        for (int k = 0; k < (int)curve2.size(); k++){
                            std::cout << "- - - - - 8" << std::endl;
                            
                            if ((j+c1+k-m) < int(cnt.size())){
                                cv::Point pkt1, pkt2;
                                pkt1 = cnt[j+c1+k-m];
                                pkt2 = cnt2[j+c1+k+m];

                                std::vector<cv::Point> line2;
                                getLinePoints(j,pkt1, pkt2, line2);
                                
                                int middle = (int)line2.size()/2;
                                std::cout << "- - - - - 9" << std::endl;
                                if (bg_img.at<uchar>(line2[middle]) == 0){
                                    
                                    cv::Point midpoint;
                                    midpoint.x = line2[middle].x;
                                    midpoint.y = line2[middle].y;
                                    
                                    double dst = sqrt(pow((curve2[k].x-midpoint.x),2.) + pow((curve2[k].y-midpoint.y),2.));
                                    
                                    if (dst  >= win_dst and dst >= min_dst){
                                        
                                        std::cout << "- - - - - 10" << std::endl;
                                        win_dst = dst;
                                        winner = cnt2[j+c1+k];
                                        winner_middle = midpoint;
                                        pos = j+c1+k;
                                    }
                                }
                            }
                        }

                        std::cout << "- - - save winning point and corresponding middle point." << std::endl;
                        
                        if (winner.x != 0 and winner.y != 0){
                            con_points.push_back(winner);
                            midpoints.push_back(winner_middle);
                            j = pos+3;
                        }
                        else j += 3;
                    }
                    else j += 3;
                }
                std::cout << "14" << std::endl;
            }

            cont_array.push_back(con_points);
            midpoint_array.push_back(midpoints);
        }
    
    }

    /**
     *  calculate all pixel points belonging to a straight line from pt1 to pt2
     * Method: rasterization of lines
     */
    void getLinePoints(int j, cv::Point pt1, cv::Point pt2, std::vector<cv::Point> &linePoints){
        
        int delta_x = pt2.x - pt1.x;	// difference x
        int delta_y = pt2.y - pt1.y;	// difference y
        double m = (double)delta_y/delta_x;		// rise

        if (abs(delta_x) >= abs(delta_y)){
            linePoints.push_back(pt1);
            cv::Point pi;
            for (int i=1; i < abs(delta_x)+1; i++){
                if (delta_x < 0){
                    pi.x = linePoints[0].x -i;
                }else{
                    pi.x = linePoints[0].x +i;
                }

                pi.y = myRounding(m*(pi.x-linePoints[0].x) + linePoints[0].y);
                linePoints.push_back(pi);
            }
        }
        else{
            //linePoints.reserve(abs(delta_y)+1);
            linePoints.push_back(pt1);
            cv::Point pi;
            for (int i=1; i < abs(delta_y)+1; i++){
                if (delta_y < 0){
                    pi.y = linePoints[0].y -i;
                }else{
                    pi.y = linePoints[0].y +i;
                }

                pi.x = singleCells::myRounding((pi.y-linePoints[0].y)/m + linePoints[0].x);
                linePoints.push_back(pi);

            }
        }
    }

    int myRounding(double val){

        int rval = (int)val;
        double diff = val-rval;

        if (diff > 0.5){
            rval += 1;
        }
        else if (diff == 0.5){
            // if odd number -> round up
            if (rval%2 != 0){
                rval += 1;
            }
        }

        return rval;
    }

    /**
     * Find the cut lines between the concavity points of one contour with the
     * help of the midpoint of the chord and the directional vector
     */
    void findCutLines(std::vector<std::vector<cv::Point>> cont_array, std::vector<std::vector<cv::Point>> midpoint_array, std::vector<std::vector<cv::Point>> contours, std::vector<cv::Vec4i> hierarchy, int min_width, int max_width, int min_length, int max_length, cv::Mat bg_img, std::vector<std::vector<std::vector<cv::Point>>> &cut_lines_array){

        for (int i = 0; i < (int)cont_array.size(); i++){
            if (hierarchy[i][3] != -1 and cont_array[i].size() > 0){
                while ((int)cont_array.size() > 0){
                    
                    std::vector<cv::Point>::iterator temp = cont_array[i].end();
                    
                    cont_array[hierarchy[i][3]].push_back(*temp);
                    cont_array[i].pop_back();

                    temp = midpoint_array[i].end();
                    midpoint_array[hierarchy[i][3]].push_back(*temp);
                    midpoint_array[i].pop_back();
                }
            }
        }

        std::vector<std::vector<cv::Point>> cut_lines;
        for (int i = 0; i < (int)cont_array.size(); i ++){
            
            cut_lines.clear();

            if (cont_array[i].size() > 1){
                for (int j = 0; j < (int)cont_array[i].size(); j++){
            
                    double width = min_width;
                    double length = min_length;
                    bool found = false;
                    
                    cv::Point cp = cont_array[i][j];
                    std::vector<cv::Point> set1 = cont_array[i];
                    set1.erase(set1.begin()+j);
                    cv::Point point;
                    
                    while (width <= max_width and found == false){
                        std::vector<cv::Point> rect;
                        getRect(cont_array[i][j], midpoint_array[i][j], width, length, rect);
                        point = findPoint(i,cp, rect, set1, found);
                        width += 2;
                    }
                    while (length <= max_length and  found == false){
                        std::vector<cv::Point> rect;
                        singleCells::getRect(cont_array[i][j], midpoint_array[i][j], width, length,rect);
                        point = singleCells::findPoint(i,cp, rect, set1, found);
                        length += 2;
                    }
                    
                    bool single = true;
                    for (int k = 0; k < (int)cut_lines.size(); k++){
                        if (point == cut_lines[k][0] or point == cut_lines[k][1]){
                            if (cp == cut_lines[k][1] or cp == cut_lines[k][0]){
                                single = false;
                            }
                        }
                    }
                    
                    bool inner = true;
                    std::vector<cv::Point> cutline;
                    singleCells::getLinePoints(j,cp, point, cutline);
                    
                    double summe = 0;
                    for (int k = 0; k < (int)cutline.size(); k++){
                        summe += bg_img.at<uchar>(cutline[k]);
                    }
                    if (summe/255 < cutline.size()-2){
                        inner = false;
                    }
                    if (single and found and inner){
                        std::vector<cv::Point> startendpoints(2);
                        
                        startendpoints[0] = cp;
                        startendpoints[1] = point;
                        cut_lines.push_back(startendpoints);
                    }
                    else if (found == false){
                        cv::Point mp = midpoint_array[i][j];
                        double width = 4;
                        double length = max_length;
                        std::vector<cv::Point> rect;
                        getRect(cp, mp, width, length, rect);

                        std::vector<double> ad, ab;
                        ad.push_back(-(rect[0].y-rect[3].y));
                        ad.push_back((rect[0].x-rect[3].x));
                        ab.push_back(-(rect[0].y-rect[1].y));
                        ab.push_back((rect[0].x-rect[1].x));

                        double l_ab = sqrt(pow(ab[0],2.)+pow(ab[1],2.));
                        double l_ad = sqrt(pow(ad[0],2.)+pow(ad[1],2.));

                        ab[0] = ab[0]/l_ab;
                        ab[1] = ab[1]/l_ab;
                        ad[0] = ad[0]/l_ad;
                        ad[1] = ad[1]/l_ad;

                        std::vector<std::vector<double>> proji;

                        for (int j=0; j < 4; j++){
                            std::vector<double> temp;
                            
                            temp.push_back((ad[0]*rect[j].x) + (ad[1]*rect[j].y));
                            temp.push_back((ab[0]*rect[j].x) + (ab[1]*rect[j].y));
                            proji.push_back(temp);
                        }

                        double min1 = std::min(proji[3][0], std::min(proji[2][0], std::min(proji[1][0], proji[0][0])));
                        double max1 = std::max(proji[3][0], std::max(proji[2][0], std::max(proji[1][0], proji[0][0])));
                        double min2 = std::min(proji[3][1], std::min(proji[2][1], std::min(proji[1][1], proji[0][1])));
                        double max2 = std::max(proji[3][1], std::max(proji[2][1], std::max(proji[1][1], proji[0][1])));

                        std::vector<cv::Point> points;

                        for (int k = 0; k < (int)contours[i].size(); k++){
                            
                            double p1 = (ad[0]*contours[i][k].x + ad[1]*contours[i][k].y);
                            double p2 = (ab[0]*contours[i][k].x + ab[1]*contours[i][k].y);
                            
                            if (min1 <= p1 and p1 <= max1 and min2 <= p2 and p2 <= max2){
                                points.push_back(contours[i][k]);
                            }
                        }
                        
                        double dst = 0;
                        cv::Point point;
                        for (int k = 0; k < (int)points.size(); k++){
                            double dst_temp = sqrt(pow((points[k].x-cp.x),2.) + pow((points[k].y-cp.y),2.));
                            if (dst_temp > dst){
                                point.x = points[k].x;
                                point.y = points[k].y;
                                dst = dst_temp;
                            }
                        }
                        
                        std::vector<cv::Point> tempv;
                        tempv.push_back(cp);
                        tempv.push_back(point);
                        cut_lines.push_back(tempv);
                    }
                }
            }
            else if (cont_array[i].size() == 1){
                cv::Point cp = cont_array[i][0];
                cv::Point mp = midpoint_array[i][0];

                double width = 4;
                double length = max_length;

                std::vector<cv::Point> rect;
                getRect(cp, mp, width, length, rect);

                std::vector<double> ad, ab;
                ad.push_back(-(rect[0].y-rect[3].y));
                ad.push_back((rect[0].x-rect[3].x));
                ab.push_back(-(rect[0].y-rect[1].y));
                ab.push_back((rect[0].x-rect[1].x));

                double l_ab = sqrt(pow(ab[0],2.)+pow(ab[1],2.));
                double l_ad = sqrt(pow(ad[0],2.)+pow(ad[1],2.));

                ab[0] = ab[0]/l_ab;
                ab[1] = ab[1]/l_ab;
                ad[0] = ad[0]/l_ad;
                ad[1] = ad[1]/l_ad;

                std::vector<std::vector<double>> proji;

                for (int j=0; j < 4; j++){
                    std::vector<double> temp;
                    
                    temp.push_back((ad[0]*rect[j].x) + (ad[1]*rect[j].y));
                    temp.push_back((ab[0]*rect[j].x) + (ab[1]*rect[j].y));
                    proji.push_back(temp);
                }

                double min1 = std::min(proji[3][0], std::min(proji[2][0], std::min(proji[1][0], proji[0][0])));
                double max1 = std::max(proji[3][0], std::max(proji[2][0], std::max(proji[1][0], proji[0][0])));
                double min2 = std::min(proji[3][1], std::min(proji[2][1], std::min(proji[1][1], proji[0][1])));
                double max2 = std::max(proji[3][1], std::max(proji[2][1], std::max(proji[1][1], proji[0][1])));

                std::vector<cv::Point> points;

                for (int k = 0; k < (int)contours[i].size(); k++){
                    double p1 = (ad[0]*contours[i][k].x + ad[1]*contours[i][k].y);
                    double p2 = (ab[0]*contours[i][k].x + ab[1]*contours[i][k].y);

                    if (min1 <= p1 and p1 <= max1 and min2 <= p2 and p2 <= max2){
                        points.push_back(contours[i][k]);
                    }
                }
                
                double dst = 0;
                cv::Point point;
                for (int k = 0; k < (int)points.size(); k++){
                    double dst_temp = sqrt(pow((points[k].x-cp.x),2.) + pow((points[k].y-cp.y),2.));
                    if (dst_temp > dst){
                        point.x = points[k].x;
                        point.y = points[k].y;
                        dst = dst_temp;
                    }
                }
                
                std::vector<cv::Point> tempv;
                tempv.push_back(cp);
                tempv.push_back(point);
                cut_lines.push_back(tempv);
            }

            cut_lines_array.push_back(cut_lines);
        
        }

    }

    /**
     * calculate the corner points a, b, c and d of a rectangle
     * use the concavity point cp, the midpoint of the chord mp and
     * length and width oft he rectangle
     */
    void getRect(cv::Point cp, cv::Point mp, double width, double length, std::vector<cv::Point> &rect){
       
        std::vector<double> a, b, c, d;
        double dst = width/2.;

        // direction vector from mp to cp
        std::vector<double> vector;			
        vector.push_back(cp.x-mp.x);
        vector.push_back(cp.y-mp.y);

        if (vector[0] == 0 and vector[1] != 0){
            
            if (vector[1] > 0){
                a.push_back(cp.x-dst);
                a.push_back(cp.y);
                b.push_back(cp.x+dst);
                b.push_back(cp.y);
                c.push_back(b[0]);
                c.push_back(b[1]+length);
                d.push_back(a[0]);
                d.push_back(a[1]+length);
            }
            else if(vector[1] < 0){
                a.push_back(cp.x+dst);
                a.push_back(cp.y);
                b.push_back(cp.x-dst);
                b.push_back(cp.y);
                c.push_back(b[0]);
                c.push_back(b[1]-length);
                d.push_back(a[0]);
                d.push_back(a[1]-length);
            }
        }
        else if (vector[1] == 0 and vector[0] != 0){
            
            if (vector[0] > 0){
                a.push_back(cp.x);
                a.push_back(cp.y+dst);
                b.push_back(cp.x);
                b.push_back(cp.y-dst);
                c.push_back(b[0]+length);
                c.push_back(b[1]);
                d.push_back(a[0]+length);
                d.push_back(a[1]);
            }
            else if(vector[0] < 0){
                a.push_back(cp.x);
                a.push_back(cp.y-dst);
                b.push_back(cp.x);
                b.push_back(cp.y+dst);
                c.push_back(b[0]-length);
                c.push_back(b[1]);
                d.push_back(a[0]-length);
                d.push_back(a[1]);
            }
        }
        else if (vector[0] != 0 and vector[1] != 0){
            
            double m1 = (double)vector[1]/(double)vector[0];
            double m2 = -1./m1;

            double var1 = sqrt((pow(dst,2.))/(1+(pow(m2,2.))));
            double var2 = sqrt((pow(length,2.))/(1+(pow(m1,2.))));

            if (vector[0] > 0){
                if (m2 < 0){
                    b.push_back(var1+cp.x);
                    b.push_back((m2*var1)+cp.y);
                    a.push_back(-var1+cp.x);
                    a.push_back(-(m2*var1)+cp.y);
                }
                else{
                    a.push_back(var1+cp.x);
                    a.push_back((m2*var1)+cp.y);
                    b.push_back(-var1+cp.x);
                    b.push_back(-(m2*var1)+cp.y);
                }
                
                d.push_back(var2+a[0]);
                d.push_back(m1*var2+a[1]);
                c.push_back(var2+b[0]);
                c.push_back(m1*var2+b[1]);
            }
            else{
                
                if (m2 < 0){
                    a.push_back(var1+cp.x);
                    a.push_back((m2*var1)+cp.y);
                    b.push_back(-var1+cp.x);
                    b.push_back(-(m2*var1)+cp.y);
                }
                else{
                    b.push_back(var1+cp.x);
                    b.push_back((m2*var1)+cp.y);
                    a.push_back(-var1+cp.x);
                    a.push_back(-(m2*var1)+cp.y);
                }

                d.push_back(-var2+a[0]);
                d.push_back(-(m1*var2)+a[1]);
                c.push_back(-var2+b[0]);
                c.push_back(-(m1*var2)+b[1]);
            }
        }

        cv::Point aa, bb, cc, dd;

        aa.x = myRounding(a[0]);
        aa.y = myRounding(a[1]);
        bb.x = myRounding(b[0]);
        bb.y = myRounding(b[1]);
        cc.x = myRounding(c[0]);
        cc.y = myRounding(c[1]);
        dd.x = myRounding(d[0]);
        dd.y = myRounding(d[1]);

        rect.push_back(aa);
        rect.push_back(bb);
        rect.push_back(cc);
        rect.push_back(dd);
    }

    /** 
     * checks, whether a point from set set1 lies in the rectangle and returns it.
     * if more than one point lies in the rectangle, the point with the smallest
     * distance to cp will be returned.
     * Method: Separating Axis Theorem (SAT)
     */
    cv::Point findPoint(int ix, cv::Point cp, std::vector<cv::Point> rect, std::vector<cv::Point> set1, bool& found ){
       
        cv::Point point;
        point.x = 0;
        point.y = 0;
        double dst = std::numeric_limits<double>::infinity();

        std::vector<double> ad, ab;

        ad.push_back(-(rect[0].y-rect[3].y));
        ad.push_back((rect[0].x-rect[3].x));
        ab.push_back(-(rect[0].y-rect[1].y));
        ab.push_back((rect[0].x-rect[1].x));

        double l_ab = sqrt(pow(ab[0],2.)+pow(ab[1],2.));
        double l_ad = sqrt(pow(ad[0],2.)+pow(ad[1],2.));

        ab[0] = ab[0]/l_ab;
        ab[1] = ab[1]/l_ab;
        ad[0] = ad[0]/l_ad;
        ad[1] = ad[1]/l_ad;

        std::vector<std::vector<double>> proji;

        for (int j=0; j < 4; j++){
            std::vector<double> temp;
            
            temp.push_back((ad[0]*rect[j].x) + (ad[1]*rect[j].y));
            temp.push_back((ab[0]*rect[j].x) + (ab[1]*rect[j].y));
            proji.push_back(temp);
        }

        double min1 = std::min(proji[3][0], std::min(proji[2][0], std::min(proji[1][0], proji[0][0])));
        double max1 = std::max(proji[3][0], std::max(proji[2][0], std::max(proji[1][0], proji[0][0])));
        double min2 = std::min(proji[3][1], std::min(proji[2][1], std::min(proji[1][1], proji[0][1])));
        double max2 = std::max(proji[3][1], std::max(proji[2][1], std::max(proji[1][1], proji[0][1])));

        for (int i = 0; i < (int)set1.size(); i++){
            double p1 = (ad[0]*set1[i].x) + (ad[1]*set1[i].y);
            double p2 = (ab[0]*set1[i].x) + (ab[1]*set1[i].y);

            double dst_temp = sqrt(pow((cp.x-set1[i].x),2.) + pow((cp.y-set1[i].y),2.));

            if (min1 <= p1 and p1 <= max1 and min2 <= p2 and p2 <= max2 and dst_temp < dst){
                found = true;
                point = set1[i];
                dst = dst_temp;
            }
        }

        return point;
    }

    /**
     * Post processing of the found lines.
     * Solving triangles, sharp angles and squares.
     */
    void postprocessCutLines(std::vector<std::vector<std::vector<cv::Point>>> &cut_lines_array){
        
        std::vector<int> degrees_num;
        std::vector<cv::Point> degrees_points;

        /* Go through all clumps separately and consider their cut-lines. */
        for (int i = 0; i < (int)cut_lines_array.size(); i++){

            /* Create and fill hash table with points from current clump according to
            * their hash value. The algorithms for the calculation of the hash values
            * and insertion of points into the table is adapted according to the
            * implementation of the python dictionary.
            * This step is necessary to order the points in the same way,
            * that they were ordered in the python version of this program .
            * This ensures, that the cut-lines are processed in the same order,
            * in which they were processed in the python program, which in turn
            * guarantees the same end results.
            *
            * Start with hash table of size 2^nn (= 8), resize when more than 2/3 of the space is filled
            * to a 4x bigger table and reallocate all previous points.
            * Initialize with points (0,0).
            * Fill hash table with points according to their hash value.
            */

            int nn = 3;
            int len_table = (int)pow(2., nn);

            /* Initialization */
            cv::Point *degrees_hash = new cv::Point[len_table];
            int filledpos = 0;
            for (int k = 0; k < len_table; k++){
                degrees_hash[k] = cv::Point(0,0);
            }

            for (int j = 0; j < (int)cut_lines_array[i].size(); j++){

                for (int k = 0; k < 2; k++){

                    cv::Point point = cut_lines_array[i][j][k];

                    long pos;
                    bool notcontained = true;
                    hashSearchPoint(degrees_hash, point, pos, notcontained, nn);

                    if (notcontained){
                        
                        if ((filledpos + 1) < ((2./3.) * len_table)){
                            degrees_hash[pos] = point;
                            filledpos += 1;
                        }
                        else{
                            cv::Point *temp_degrees_hash = new cv::Point[len_table*4];
                            nn += 2;
                            filledpos = 0;
                        
                            for (int k = 0; k < len_table*4; k++){
                                temp_degrees_hash[k] = cv::Point(0,0);
                            }
                        
                            for (int k = 0; k < len_table; k++){
                                if (degrees_hash[k].x != 0 and degrees_hash[k].y != 0){
                                    cv::Point temp_point = degrees_hash[k];
                                    hashSearchPoint(temp_degrees_hash, temp_point, pos, notcontained, nn);
                                    temp_degrees_hash[pos] = temp_point;
                                    filledpos += 1;
                                }
                            }
                            len_table = len_table*4;

                            delete [] degrees_hash;
                            degrees_hash = new cv::Point[len_table];
                            degrees_hash = temp_degrees_hash;
                            hashSearchPoint(degrees_hash, point, pos, notcontained, nn);
                            degrees_hash[pos] = point;
                            filledpos += 1;
                        }
                    }
                }
            }

            degrees_num.clear();
            degrees_points.clear();

            /// Go through points and save them in degrees_points.

            for (int ii = 0; ii < (pow(2.,nn)); ii++){
                if (degrees_hash[ii].x != 0 and degrees_hash[ii].y != 0){				// wenn Eintrag an dieser Stelle
                    degrees_points.push_back(degrees_hash[ii]);
                }
            }

            for (int j = 0; j < (int)degrees_points.size(); j++){
                cv::Point temppoint = degrees_points[j];
                int count = 0;

                for(int k = 0; k < (int)cut_lines_array[i].size(); k++){
                    for (int ii = 0; ii < 2; ii++){
                        if(cut_lines_array[i][k][ii].x == temppoint.x and cut_lines_array[i][k][ii].y == temppoint.y){
                            count += 1;
                        }
                    }
                }

                degrees_num.push_back(count);
            }

            std::vector<int> degree2 = degrees_num;

            /* sort the point vector according to the highest degree
            * calculate sum over all degrees_num
            */

            std::vector<cv::Point> points_sorted;
            int summe = 0;
            for (int l = 0; l < (int)degree2.size(); l++){
                summe += degree2[l];
            }

            while (summe > 0){
                std::vector<int>::iterator max = degree2.begin();
                int max_pos = 0;
                int pos = 0;
                
                for (std::vector<int>::iterator it = degree2.begin(); it != degree2.end(); it++, pos++){
                    if (*it > *max){
                        max = it;
                        max_pos = pos;
                    }
                }
                
                points_sorted.push_back(degrees_points[max_pos]);
                summe = summe - *max;
                degree2[max_pos] = 0;

            }

            /// go through all the points and check, if lines should be deleted

            for(int j = 0; j < (int)points_sorted.size(); j++){
               
                // the current point
                cv::Point cp = points_sorted[j];								
                bool triangle = false;
                // positions of relevant lines in the line vector
                std::vector<int> pos;	
                // points connected to the current point											
                std::vector<cv::Point> points;										

                int it, it2, it3;
                // get degree of current point
                int iter = -1;														
                for (int l = 0; l < (int)degrees_points.size(); l++){
                    if (degrees_points[l].x == cp.x and degrees_points[l].y == cp.y){
                        iter = l;
                    }
                }

                if (iter != -1 and degrees_num[iter] == 2){
                    // get connected points
                    points = getPoints(cp, cut_lines_array[i], triangle, pos);	
                    cv::Point c;
                    c.x = myRounding((double)1/3*(cp.x+points[0].x+points[1].x));
                    c.y = myRounding((double)1/3*(cp.y+points[0].y+points[1].y));
                    if (triangle){

                        // consider all related points and set the centre c
                        for (int k = 0; k < (int)points.size(); k++){			
                            cut_lines_array[i][pos[k]][0] = points[k];
                            cut_lines_array[i][pos[k]][1] = c;
                            
                            // update degree of the points 
                            it = -1;											
                            for (int l = 0; l < (int)degrees_points.size(); l++){
                                if (degrees_points[l].x == points[k].x and degrees_points[l].y == points[k].y){
                                    it = l;
                                }
                            }
                            if (it != -1){
                                degrees_num[it] -= 1;
                            }
                        }

                        // also update the current point and its corresponding degree
                        cut_lines_array[i][pos[2]][0] = cp;					
                        cut_lines_array[i][pos[2]][1] = c;

                        // get degree of current point
                        it = -1;														
                        for (int l = 0; l < (int)degrees_points.size(); l++){
                            if (degrees_points[l].x == cp.x and degrees_points[l].y == cp.y){
                                it = l;
                            }
                        }
                        if (it != -1){
                            degrees_num[it] -= 1;
                        }

                        // save degree of c 
                        degrees_points.push_back(c);				
                        degrees_num.push_back(3);
                    }
                    else{
                        
                        for (int k = 0; k< (int)points.size(); k++){
                            cut_lines_array[i][pos[k]][0] = points[k];
                            cut_lines_array[i][pos[k]][1] = c;
                        }
                        
                        // add new line 
                        std::vector<cv::Point> new_cut_line;							
                        new_cut_line.push_back(cp);
                        new_cut_line.push_back(c);
                        cut_lines_array[i].push_back(new_cut_line);

                        // get degree of current point
                        it = -1;														
                        for (int l = 0; l < (int)degrees_points.size(); l++){
                            if (degrees_points[l].x == cp.x and degrees_points[l].y == cp.y){
                                it = l;
                            }
                        }
                        if (it != -1){
                            degrees_num[it] -= 1;
                        }
                        degrees_points.push_back(c);
                        degrees_num.push_back(3);
                    }
                }
                else if(iter != -1 && degrees_num[iter] == 3){

                    // get all connected points
                    points = getPoints(cp, cut_lines_array[i], triangle, pos);	
                    
                    std::vector<cv::Point> vec;
                    // direction vectors from cp to all points
                    for(int k = 0; k < (int)points.size(); k++){				
                        cv::Point temp;
                        temp.x = points[k].x-cp.x;
                        temp.y = points[k].y -cp.y;
                        vec.push_back(temp);
                    }
                    // first direction vector again at the end to get all neighboring lines
                    vec.push_back(vec[0]);										

                    // angle vector
                    std::vector<double> angles;								
                    // calculate the angles between neighboring lines		
                    for (int k = 1; k < (int)vec.size(); k++){				
                        double vecprod = vec[k-1].x*vec[k].x + vec[k-1].y*vec[k].y;
                        double l1 = sqrt(pow(vec[k].x,2.) + pow(vec[k].y,2.));
                        double l2 = sqrt(pow(vec[k-1].x,2.) + pow(vec[k-1].y,2.));
                        angles.push_back(acos(vecprod/(l1*l2))*180/M_PI);
                    }

                    int p = 0;
                    // max angle
                    double ang = 360;	

                    // find the smallest angle and save its index in the angle vector
                    for (int k = 0; k < (int)angles.size(); k++){				
                        if (angles[k] < ang){
                            ang = angles[k];
                            p = k;
                        }
                    }

                    cv::Point p1, p2;
                    // get the end points for the two lines enclosing the smallest angle
                    if (p == (int)points.size()-1){								
                        p1 = points[p];
                        p2 = points[0];
                    }
                    else{
                        p1 = points[p];
                        p2 = points[p+1];
                    }

                    std::vector<std::vector<cv::Point>>::iterator iterat;
                    // get the degree of these two points
                    it = -1;															
                    it2 = -1;
                    // get degree of current point
                    int it = -1;														
                    for (int l = 0; l < (int)degrees_points.size(); l++){
                        if (degrees_points[l].x == p1.x and degrees_points[l].y == p1.y){
                            it = l;
                        }
                    }
                    for (int l = 0; l < (int)degrees_points.size(); l++){
                        if (degrees_points[l].x == p2.x and degrees_points[l].y == p2.y){
                            it2 = l;
                        }
                    }

                    // remove one of the lines (biggest angle)
                    if (it != -1 and degrees_num[it] == 2 and it2 != -1 and degrees_num[it2] == 2){		
                        std::vector<cv::Point> templine;
                        if (p-1 < 0){
                            templine.push_back(points[points.size()-1]);
                        }
                        else {
                            templine.push_back(points[p-1]);
                        }                         
                        
                        templine.push_back(cp);

                        for (int ij = 0; ij < (int)cut_lines_array[i].size(); ij++){
                            if((cut_lines_array[i][ij][0].x == templine[0].x and cut_lines_array[i][ij][0].y == templine[0].y and cut_lines_array[i][ij][1].x == templine[1].x and cut_lines_array[i][ij][1].y == templine[1].y) or
                                    (cut_lines_array[i][ij][0].x == templine[1].x and cut_lines_array[i][ij][0].y == templine[1].y and cut_lines_array[i][ij][1].x == templine[0].x and cut_lines_array[i][ij][1].y == templine[0].y)){
                                cut_lines_array[i].erase(cut_lines_array[i].begin()+ij);
                            }
                        }

                        // update degrees of the point and cp
                        it3 = -1;
                        for (int l = 0; l < (int)degrees_points.size(); l++){
                            if (degrees_points[l].x == templine[0].x and degrees_points[l].y == templine[0].y){
                                it3 = l;
                            }
                        }
                        if (it3 != -1){
                            degrees_num[it3] -= 1;
                        }
                        it3 = -1;
                        for (int l = 0; l < (int)degrees_points.size(); l++){
                            if (degrees_points[l].x == templine[1].x and degrees_points[l].y == templine[1].y){
                                it3 = l;
                            }
                        }
                        if (it3 != -1){
                            degrees_num[it3] -= 1;
                        }
                    }
                    else if (it != -1 and degrees_num[it] == 2 and it2 != -1 and degrees_num[it2] < 2){
                        // remove line cp-p1
                        std::vector<cv::Point> templine;
                        templine.push_back(cp);
                        templine.push_back(p1);
                        
                        for (int ij = 0; ij < (int)cut_lines_array[i].size(); ij++){
                            if((cut_lines_array[i][ij][0].x == cp.x and cut_lines_array[i][ij][0].y == cp.y and cut_lines_array[i][ij][1].x == p1.x and cut_lines_array[i][ij][1].y == p1.y) or
                                    (cut_lines_array[i][ij][0].x == p1.x and cut_lines_array[i][ij][0].y == p1.y and cut_lines_array[i][ij][1].x == cp.x and cut_lines_array[i][ij][1].y == cp.y)){
                                cut_lines_array[i].erase(cut_lines_array[i].begin()+ij);
                            }
                        }

                        // reduce the degree of cp and p1 by one
                        it3 = -1;
                        for (int l = 0; l < (int)degrees_points.size(); l++){
                            if (degrees_points[l].x == templine[0].x and degrees_points[l].y == templine[0].y){
                                it3 = l;
                            }
                        }
                        if (it3 != -1){
                            degrees_num[it3] -= 1;
                        }
                        it3 = -1;
                        for (int l = 0; l < (int)degrees_points.size(); l++){
                            if (degrees_points[l].x == templine[1].x and degrees_points[l].y == templine[1].y){
                                it3 = l;
                            }
                        }
                        if (it3 != -1){
                            degrees_num[it3] -= 1;
                        }

                    }
                    else if (it != -1 and degrees_num[it] < 2 and it2 != -1 and degrees_num[it2] == 2){
                        
                        // remove line cp-p2
                        std::vector<cv::Point> templine;
                        templine.push_back(cp);
                        templine.push_back(p2);

                        for (int ij = 0; ij < (int)cut_lines_array[i].size(); ij++){
                            if((cut_lines_array[i][ij][0].x == cp.x and cut_lines_array[i][ij][0].y == cp.y and cut_lines_array[i][ij][1].x == p2.x and cut_lines_array[i][ij][1].y == p2.y) or
                                    (cut_lines_array[i][ij][0].x == p2.x and cut_lines_array[i][ij][0].y == p2.y and cut_lines_array[i][ij][1].x == cp.x and cut_lines_array[i][ij][1].y == cp.y)){
                                cut_lines_array[i].erase(cut_lines_array[i].begin()+ij);
                            }
                        }

                        // reduce the degree of cp and p1 by one
                        it3 = -1;
                        for (int l = 0; l < (int)degrees_points.size(); l++){
                            if (degrees_points[l].x == templine[0].x and degrees_points[l].y == templine[0].y){
                                it3 = l;
                            }
                        }
                        if (it3 != -1){
                            degrees_num[it3] -= 1;
                        }
                        it3 = -1;
                        for (int l = 0; l < (int)degrees_points.size(); l++){
                            if (degrees_points[l].x == templine[1].x and degrees_points[l].y == templine[1].y){
                                it3 = l;
                            }
                        }
                        if (it3 != -1){
                            degrees_num[it3] -= 1;
                        }

                    }
                    else if (it != -1 and degrees_num[it] < 2 and it2 != -1 and degrees_num[it2] < 2){
                        /* do normal postprocessing between four points (rectangle instead of triangle)
                        * calculate middle point c
                        */
                        cv::Point c;
                        c.x = myRounding((double)1/4*(cp.x+points[0].x+points[1].x+points[2].x));
                        c.y = myRounding((double)1/4*(cp.y+points[0].y+points[1].y+points[2].y));
                        // replace cp by c in all three lines

                        for (int k = 0; k < (int)points.size(); k++){
                            cut_lines_array[i][pos[k]][0] = points[k];
                            cut_lines_array[i][pos[k]][1] = c;
                        }
                        
                        // create new line between cp and c
                        std::vector<cv::Point> new_line;
                        new_line.push_back(cp);
                        new_line.push_back(c);
                        cut_lines_array[i].push_back(new_line);
                        
                        // update degrees of cp and c
                        it3 = -1;
                        for (int l = 0; l < (int)degrees_points.size(); l++){
                            if (degrees_points[l].x == cp.x and degrees_points[l].y == cp.y){
                                it3 = l;
                            }
                        }

                        if (it3 != -1){
                            degrees_num[it3] -= 2;
                        }

                        degrees_points.push_back(c);
                        degrees_num.push_back(4);

                    }
                }
            }
        }
    }

    /**
     * find all points that are connected to the current point via a line from lines.
     * check whether they build a tringle.
     * save position of the line in the array for later?
     */
    std::vector<cv::Point> getPoints(cv::Point point, std::vector<std::vector<cv::Point>> lines, bool& triangle, std::vector<int>& num){
        
        std::vector<cv::Point> points;
        // count variable for points and num
        int n = 0;									
        for (int i= 0; i <(int)lines.size(); i++){
            
            if (lines[i][0].x == point.x and lines[i][0].y == point.y){
                points.push_back(lines[i][1]);
                num.push_back(i);
                n += 1;
            }
            else if (lines[i][1].x == point.x and lines[i][1].y == point.y){
                points.push_back(lines[i][0]);
                num.push_back(i);
                n += 1;
            }

        }

        for (int i = 0; i < (int)lines.size(); i++){
            
            if (((points[0].x == lines[i][0].x and points[0].y == lines[i][0].y)
                    or (points[0].x == lines[i][1].x and points[0].y == lines[i][1].y))
                    and ((points[1].x == lines[i][0].x and points[1].y == lines[i][0].y)
                            or (points[1].x == lines[i][1].x and points[1].y == lines[i][1].y))){
                triangle = true;
                num.push_back(i);
                n += 1;
            }

        }

        return points;
    }

    /**
     *  calculate hash value from point and convert to binary 
     */
    void hashSearchPoint(cv::Point  hashtable[], cv::Point point, long &pos, bool &notcontained, int &nn){

        long hash_val = hashvalue(point);
        std::vector<int> binary;
        binary = decimal2binary(hash_val);

        // calculate position of point in hash table from the last nn positions of the binary number.
        
        reverse(binary.begin(), binary.end());
        pos = 0;
        notcontained = true;

        for (int it=0; it < nn; it++){
            pos += binary[it]*(pow(2.,it));
        }

        /* If there is no other point at position pos, return this position.
        * Otherwise test, if position already contains the same point or
        * if it contains another. If it contains the same point, set notcontained
        * to FALSE. If it contains another point, find the next free spot where the
        * point can be inserted. Do this until you find a free spot or the point
        * itself. Calculation of the next spot according to implementation for the
        * python dictionary.
        */

        if (hashtable[pos].x != 0 or hashtable[pos].y != 0){
        
            if (hashtable[pos].x != point.x or hashtable[pos].y != point.y){
         
                bool notfound = true;
                long perturb = hash_val;
         
                while (notfound) {
                    pos = (5*pos)+1+perturb;
                    perturb >>= 5;
                    pos = pos % (int)(pow(2.,nn));
                    
                    if (hashtable[pos].x == point.x and hashtable[pos].y == point.y){
                        notfound = false;
                        notcontained = false;
                    }
                    else if (hashtable[pos].x == 0 and hashtable[pos].y == 0){
                        notfound = false;
                        notcontained = true;
                    }
                }
            }
            else {
                notcontained = false;
            }
        }
        else {
            notcontained = true;
        }

    }

    /**
     * Calculates a hashvalue from a point.
	 * Uses the same calculations that the builtin hash function in python 2.7.3
	 * uses for the calculation of a hash value from a tuple.
	 * See method tuplehash() in file tupleobject.c in the python sources.
	 */
    long hashvalue(cv::Point punkt){
	
        // Initialization
        long hashwert = 0x345678L;
        long mult = 1000003L;
        int length = 1;

        // Calculation with x and y values of the point. Order matters.
        hashwert = (hashwert ^ punkt.x) * mult;
        mult += (long)(82520L + length + length);
        length--;
        hashwert = (hashwert ^ punkt.y) * mult;
        // the cast might truncate length; that doesn't change the hash stability
        mult += (long)(82520L + length + length);

        hashwert += 97531L;

        return hashwert;

    }

    std::vector<int> decimal2binary(long int z){

        std::vector<int> bitstring;

        while (z != 0){
            bitstring.push_back(z%2);
            z = z/2;
        }

        reverse(bitstring.begin(), bitstring.end());

        return bitstring;

    }

    void singleCells_segmentation(const std::vector<cv::Mat> &images, const int &N_THREADS, const std::string &OUT_BINARY, const bool &FLAG_DEBUG){

        // multi-threaded processing of images
        #pragma omp parallel for firstprivate(images) num_threads(N_THREADS)
        for(size_t i = 0; i < images.size(); ++i)
        {            
            cv::Mat binary = singleCells::processSingle(images[i], FLAG_DEBUG);

            // save binary image including PMNs
            io::write_image(OUT_BINARY, binary, i, true);

        }
    


    }

} // singleCells