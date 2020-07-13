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

#include <opencv2/opencv.hpp>


/* Function Headers in namespace "visual" */
namespace visual {
    
    /* enum for specify image operation */
    enum imgOperation {
        globThreshold = 0,
        morphology = 1,
        adaptiveThreshold = 2,
        canny = 3,
        houghTransformation = 4,
        findContours = 5,
        medianFiler = 6,
        distanceTransform = 7,
        moments = 8
    };

    void createGUI( cv::Mat &img, const imgOperation &op );

}
