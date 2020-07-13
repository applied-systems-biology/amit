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

#include "RegionSplitting.h"
#include <lemon/list_graph.h>
#include <lemon/matching.h>


namespace RegionSplitting
{   

    /**
     *
     * Region cluster: Region that will be separated
     * vector<Region> subregions: overlapping regions from previous timepoint -> used as seed
     * DO NOT CHANGE REGION IDs!
     */
    void separateClusterGMM(Region &cluster, std::vector<Region> &subregions){

        int K = subregions.size();

        cv::Mat image = cluster.getImage();
        cv::Mat tmp;
        image.copyTo(tmp);

        int minx, miny, maxx, maxy;
        cluster.getMinMax(minx, miny, maxx, maxy);

        std::vector<int> xvalues;
        std::vector<int> yvalues;

        // obersvations: x- and y-coordinates
        for(size_t i = 0; i < image.cols; i++){
            for(size_t j = 0; j < image.rows; j++){
                if(image.at<uchar>(j,i) > 0){
                    xvalues.push_back(i);
                    yvalues.push_back(j);
                }
            }
        }

        int size = xvalues.size();
        // commented out by Philipp: no usage of cv::Point center(x,y)
//        int x = 0, y = 0;
//        for(size_t i = 0; i < size; i++){
//            x += xvalues.at(i);
//            y += yvalues.at(i);
//        }
//        x /= size;
//        y /= size;
//
//        cv::Point center (y,x);

        // declare EM object and initialize model parameters
        cv::Ptr<cv::ml::EM> em_model = cv::ml::EM::create();
        em_model->setClustersNumber(K); 
        em_model->setCovarianceMatrixType(cv::ml::EM::COV_MAT_GENERIC);
        em_model->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER+cv::TermCriteria::EPS, 500, 0.1));

        // two-dimensional samples
        cv::Mat samples = cv::Mat( size, 2, CV_32FC1); 
        
        for(size_t i = 0; i < size; i++){
            samples.at<float>(i,0) = (float) xvalues.at(i);
            samples.at<float>(i,1) = (float) yvalues.at(i);
        }

        cv::Mat labels;
        cv::Mat probs;

        em_model->trainEM( samples, cv::noArray(), labels, probs );

        cv::Mat means = em_model->getMeans();
    
        std::vector<Region> tmpregions(K);

        for(size_t i = 0; i < size; i++){
            cv::Point p(xvalues.at(i), yvalues.at(i));
            int l = labels.at<int>(i);

            cv::Vec2i v(yvalues.at(i), xvalues.at(i));
            tmpregions.at(l).region_pixels.push_back(v);
            tmp.at<uchar>(p) = 255 / K*(l+1);
        }

        for(int i = 0; i < K; i++){
            cv::Point c1 (cv::saturate_cast<int>(means.at<double>(i,0)), cv::saturate_cast<int>(means.at<double>(i,1)));
            cv::circle(tmp, c1, 2, cv::Scalar(0),2);
            tmpregions.at(i).computeCentroid();
        }

        std::vector<Region>::iterator sub_it;
        for(sub_it = tmpregions.begin(); sub_it != tmpregions.end(); sub_it++){
            
            std::vector<cv::Vec2i>::iterator regpix_it;
            
            for(regpix_it = sub_it->region_pixels.begin(); regpix_it != sub_it->region_pixels.end(); regpix_it++){
                regpix_it->val[0] += minx;
                regpix_it->val[1] += miny;
            }

            sub_it->computeContour();
            sub_it->computeCentroid();
        }

//        cv::imshow("tmp", tmp);
//        cv::waitKey(0);
        
        RegionSplitting::associateCellsBipGraphMatching(K, subregions, tmpregions);

        subregions = tmpregions;
    }

    void separateClusterGMM2(RegionP &cluster, std::vector<RegionP> &subregions){

        int K = subregions.size();
    
        cv::Mat image = cluster.getImage();
        cv::Mat tmp;
        image.copyTo(tmp);

        int minx, miny, maxx, maxy;
        cluster.getMinMax(minx, miny, maxx, maxy);

        std::vector<int> xvalues;
        std::vector<int> yvalues;

        /// obersvations: x- and y-coordinates
        for(int i = 0; i < image.cols; i++){
            for(int j = 0; j < image.rows; j++){
                if(image.at<uchar>(j,i) > 0){
                    xvalues.push_back(i);
                    yvalues.push_back(j);
                }
            }
        }

        int size = xvalues.size();
        int x = 0, y = 0;
        for(int i = 0; i < size; i++){
            x += xvalues.at(i);
            y += yvalues.at(i);
        }
        x /= size;
        y /= size;

        cv::Point center (y,x);

        // declare EM object and initialize model parameters
        cv::Ptr<cv::ml::EM> em_model = cv::ml::EM::create();
        em_model->setClustersNumber(K); 
        em_model->setCovarianceMatrixType(cv::ml::EM::COV_MAT_GENERIC);
        em_model->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER+cv::TermCriteria::EPS, 500, 0.1));
        
        /// two-dimensional samples
        cv::Mat samples( size, 2, CV_32FC1); 

        for(int i = 0; i < size; i++){
            samples.at<float>(i,0) = (float) xvalues.at(i);
            samples.at<float>(i,1) = (float) yvalues.at(i);
        }

        cv::Mat labels;

        em_model->trainEM( samples, cv::noArray(), labels, cv::noArray() );

        cv::Mat means = em_model->getMeans();

        std::vector<RegionP> tmpregions(K);

        for(int i = 0; i < size; i++){
            cv::Point p(xvalues.at(i), yvalues.at(i));
            int l = labels.at<int>(i);

            cv::Vec2i v(yvalues.at(i), xvalues.at(i));
            tmpregions.at(l).region_pixels.push_back(v);
            tmp.at<uchar>(p) = 255/K*(l+1);
        }

        for(int i = 0; i < K; i++){
            cv::Point c1 (cv::saturate_cast<int>(means.at<double>(i,0)), cv::saturate_cast<int>(means.at<double>(i,1)));
            cv::circle(tmp, c1, 2, cv::Scalar(0),2);
            tmpregions.at(i).computeCentroid();
        }

        std::vector<RegionP>::iterator sub_it;
        for(sub_it = tmpregions.begin(); sub_it != tmpregions.end(); sub_it++){
            
            std::vector<cv::Vec2i>::iterator regpix_it;
            for(regpix_it = sub_it->region_pixels.begin(); regpix_it != sub_it->region_pixels.end(); regpix_it++){
                regpix_it->val[0] += minx;
                regpix_it->val[1] += miny;
            }
            sub_it->computeContour();
            sub_it->computeCentroid();
        }

        /// spread green region pixels
        std::vector<cv::Vec2i>::iterator greenit, regit;
        for(greenit = cluster.green_region_pixels.begin(); greenit != cluster.green_region_pixels.end(); greenit++){
            for(sub_it = tmpregions.begin(); sub_it != tmpregions.end(); sub_it++){
                for(regit = sub_it->region_pixels.begin(); regit != sub_it->region_pixels.end(); regit++){

                    if(greenit->val[0] == regit->val[0] && greenit->val[1] == regit->val[1]){
                        sub_it->addGreenPx(*greenit);
                        break;
                    }
                }
            }
        }

        RegionSplitting::associateCellsBipGraphMatching2(K, subregions, tmpregions);

        subregions = tmpregions;

    }

    void associateCellsBipGraphMatching(const int &K, std::vector<Region> &regions, std::vector<Region> &newregions){

        // generate cost function (pairwise Euklidean distance)
        cv::Mat cost;
        std::vector<Region*> regions_p, newregions_p;
        vecToVecP(regions, regions_p);
        vecToVecP(newregions, newregions_p);

        RegionSplitting::generateCostFunctionPairwiseEuklideanDistance(cost, regions_p, newregions_p);
    
        // associate correct cells -> use weigthed matching of bipartite graphs for solving
        std::vector<int> matching;
        RegionSplitting::minimumWeightedGraphMatching(cost, matching);

        std::vector<Region> tmp(K);
        for(size_t i = 0; i < K ; i++){
            int index = matching.at(i);
    
            tmp.at(index) = newregions.at(i);
            tmp.at(index).computeContour();
            tmp.at(index).computeCentroid();
        }

        for(size_t i = 0; i < K ; i++){
            tmp.at(i).setId(regions.at(i).getId());
        }

        newregions = tmp;
    }

    void associateCellsBipGraphMatching2(int K, std::vector<RegionP> &regions, std::vector<RegionP> &newregions){

        bool green =false;
        for(int i = 0; i < K ; i++){
            if(regions.at(i).getFungalRatioGreen() > 0){
                green = true;
                break;
            }
        }

        /// generate cost function (pairwise Euklidean distance)
        cv::Mat cost;

        if(green){
            RegionSplitting::generateCostFunctionPairwiseEuklideanDistanceAndPhagocytosis(cost, regions, newregions);
        }
        else{
            std::vector<Region*> regions_p, newregions_p;
            
            for(size_t i = 0; i < regions.size(); i++){
                Region * r = & regions.at(i);
                regions_p.push_back(r);
                
                Region * p = & newregions.at(i);
                newregions_p.push_back(p);
            }
            
            RegionSplitting::generateCostFunctionPairwiseEuklideanDistance(cost, regions_p, newregions_p);
        }

        /// associate correct cells -> use weigthed matching of bipartite graphs for solving
        std::vector<int> matching;
        RegionSplitting::minimumWeightedGraphMatching(cost, matching);

        std::vector<RegionP> tmp(K);
        for(int i = 0; i < K ; i++){
            int index = matching.at(i);
    
            tmp.at(index) = newregions.at(i);
            tmp.at(index).computeContour();
            tmp.at(index).computeCentroid();
        }

        for(int i = 0; i < K ; i++){
            tmp.at(i).setId(regions.at(i).getId());
        }

        newregions = tmp;
    }

    void generateCostFunctionPairwiseEuklideanDistance(cv::Mat &cost, std::vector<Region*> & regions1, std::vector<Region*> & regions2){

        if(regions1.size() != regions2.size())
            std::cout << "Error in generateCostFunctionPairwiseEuklideanDistance(): input region vectors do not have the same size! " << std::endl;

        int K = regions1.size();

        cv::Mat c(K, K, CV_32FC1, cv::Scalar(0));

        std::vector<Region*>::iterator regit1, regit2;

        int i = 0;
        for(regit1 = regions1.begin(); regit1 != regions1.end(); regit1++){
            int j = 0;
            for(regit2 = regions2.begin(); regit2 != regions2.end(); regit2++){
                c.at<float>(i,j) = (float) computeDistance(*regit1, *regit2);
                j++;
            }
            i++;
        }

        c.copyTo(cost);

    }

    void generateCostFunctionPairwiseEuklideanDistanceAndPhagocytosis(cv::Mat &cost, std::vector<RegionP> &regions1, std::vector<RegionP> &regions2){
	    
        if(regions1.size() != regions2.size()){
            std::cout << "Error in generateCostFunctionPairwiseEuklideanDistanceAndPhagocytosis(): input region vectors do not have the same size! " << std::endl;
        }

        int K = regions1.size();

        cv::Mat c(K, K, CV_32FC1, cv::Scalar(0));
        cv::Mat frs(K, K, CV_32FC1, cv::Scalar(0));

        std::vector<RegionP>::iterator regit1, regit2;

        int i = 0;
        for(regit1 = regions1.begin(); regit1 != regions1.end(); regit1++){
            
            int j = 0;
            
            double fr1 = regit1->getFungalRatioGreen();
            for(regit2 = regions2.begin(); regit2 != regions2.end(); regit2++){
                double fr2 = regit2->getFungalRatioGreen();
                double fr = std::abs(fr1-fr2);
                if(fr < 0.01) {
                    fr = 0.01;
                }
    
                frs.at<float>(i,j) = (float) fr;
                c.at<float>(i,j) = (float) computeDistance(*regit1, *regit2);
                j++;
            }
            i++;
        }

        cv::multiply(frs, c, c);

        c.copyTo(cost);

    }

    void minimumWeightedGraphMatching(cv::Mat &cost, std::vector<int> &match){

        match.resize(0);

        if(cost.cols != cost.rows){
            std::cout << "Error! at minimumWeightedGraphMatching(): distance matrix is not quadratic! " << std::endl;
        }
        else{
            int K = cost.cols;

            // declare graph object
            lemon::ListGraph g; 

            // add 2K nodes to the graph
            for(size_t i = 0; i < 2*K ; i++){
                g.addNode();
            }

            // add K*K edges to the graph
            for(size_t i = 0; i < K; i++){
                for(size_t j = K; j < 2*K; j++){
                    g.addEdge(g.nodeFromId(i), g.nodeFromId(j));
                }
            }

            double min, max;
            cv::minMaxLoc(cost, &min, &max);

            // add K*K Weights to an EdgeMap (type: float)
            lemon::ListGraph::EdgeMap<float> weights(g);
            int index = 0;
            for(size_t i = 0; i < K ; i++){
                for(size_t j = 0; j < K ; j++){
                    // take negative weights to determine the minimum weight matching
                    weights[g.edgeFromId(index)] = (float) max + 1.f - cost.at<float>(j,i); 
    				
                    index++;
                }
            }

            // create object for maximum weighted graph matching
            lemon::MaxWeightedMatching<lemon::ListGraph, lemon::ListGraph::EdgeMap<float> > mwm(g, weights);
            // run the algorithm
            mwm.run(); 

            // return the matchings
            int k = 0;
            for(size_t i = 0; i < K*K; i++){
                lemon::ListGraph::Edge e = g.edgeFromId(i);
                
                if(mwm.matching(e)){
                    match.push_back(i % K);
                }
                
                k++;
                
                if(k == K){
                    k = 0;
                }
            }
        }

    }

} // RegionSplitting