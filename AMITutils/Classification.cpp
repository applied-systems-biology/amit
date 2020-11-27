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

#include "Classification.h"
#include <iostream>
#include <fstream>
#include <assert.h>
#include "IOputs.h"


namespace Classification 
{

    void trainSVMonIBPdetection(const std::vector<float> &areas, const std::vector<int> &labels, cv::Ptr<cv::ml::SVM> &svm) {

        std::cout << "\tClassification::trainSVMonIBPdetection:: data - label size: " << areas.size() << " - " << labels.size() << std::endl; 

        /// number of samples has to be equal to number of labels
        assert(areas.size() == labels.size());

        // set up training data
        float trainingData[ areas.size() ];
        std::copy(areas.begin(), areas.end(), trainingData);
        
        int labels_test[ labels.size() ];
        std::copy(labels.begin(), labels.end(), labels_test);

        /// copy to cv::Mat data type
        cv::Mat trainingDataMat( areas.size() , 1, CV_32F, trainingData);
        cv::Mat labelsMat( labels.size() , 1, CV_32SC1, labels_test);

        // train the SVM
        svm = cv::ml::SVM::create();
        svm->setType(cv::ml::SVM::C_SVC);
        svm->setKernel(cv::ml::SVM::LINEAR);
        svm->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER, 300, 1e-6));
        svm->train(trainingDataMat, cv::ml::SampleTypes::ROW_SAMPLE, labelsMat);
        
        // std::cout << "SVM - SupportVectors: " << svm->getSupportVectors() << std::endl; 

        /// predict with SVM on 5 exemplary area sizes and print it
        std::vector<float> areas_test{ 2500. , 3800. , 4200. , 4600. , 5000. };

        for (size_t i = 0; i < areas_test.size(); i++) {

            float testData = areas_test[i];
            cv::Mat testDataMat( 1, 1, CV_32F, testData);

            float prediction = svm->predict( testDataMat );

            std::cout << "SVM prediction for size:\t" << areas_test[i] << " -> " << prediction << std::endl;
        } 

    }


    /**
     * This function approximates a Gaussian Mixture Distribution (GMM) to the area of a set of regions.
     *
     *
     * The parameters of the GMM are estimated with EM using the opencv library em (https://docs.opencv.org/3.0-beta/modules/ml/doc/expectation_maximization.html).
     * Labels are returned for classes in increasing order (first class has the smalles mean area, last class has the biggest mean area).
     * Labeling results is saved in the output directory path.
     *
     * @param K number of Gaussian functions
     * @param regions two-dimensional set of regions
     * @param labels class labels
     * @param outdir output directory path
     *
     */
    void approxMixDistrEM(const int &n_classes, std::vector<std::vector<Region>> &regions, cv::Mat &labels, cv::Ptr<cv::ml::EM> &em_model, const std::string &outdir, bool isFlat){

        // fill samplevector with observated areas
        std::vector<int> samplevector;
        for(size_t i = 0; i < regions.size(); i++){
            for(size_t j = 0; j < regions.at(i).size(); j++){
                int  a = regions.at(i).at(j).getArea();
                samplevector.push_back(a);
            }
        }

        int nsamples = samplevector.size();
        // one-dimensional samples
        cv::Mat samples( nsamples, 1, CV_32FC1); 

        for(size_t i = 0; i < nsamples; i++){
            samples.at<float>(i) = (float) samplevector.at(i);
        }

        /// if em_model is empty, create a new object
        if ( em_model.empty() ){
            em_model = cv::ml::EM::create();
        }

        // declare EM object and initialize parameter for EM object
        em_model->setClustersNumber(n_classes); 
        em_model->setCovarianceMatrixType(cv::ml::EM::COV_MAT_SPHERICAL);
        em_model->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER+cv::TermCriteria::EPS, 500, 0.1));
        
        // output of EM-model
        cv::Mat probs;
        cv::Mat likelihood;

        // learn gmm-model
        em_model->trainEM( samples, likelihood, labels, probs );

        // determine if labels are in the right order (0,1,2,...)
        cv::Mat means = em_model->getMeans();
        cv::Mat means_sort;
        cv::sort(means, means_sort, cv::SORT_EVERY_COLUMN + cv::SORT_ASCENDING );

        // check if classes are in an ascending order
        cv::Mat means_sub;
        cv::subtract(means, means_sort, means_sub);
        means_sub = cv::abs(means_sub);

        double sum = 0;
        for(size_t i = 0; i < n_classes; i++){
            sum += means_sub.at<double>(i,0); 
        }       

        // change labeling if it is not in an increasing order
        if(sum > 0.1) {
            // lookuptable of indizes
            std::vector<int> indizes(n_classes,0); 

            // fill lookup-table
            for(size_t i = 0; i < n_classes; i++){
                double temp = means.at<double>(i);
                for(size_t j = 0; j < n_classes; j++){
                    if(abs(temp - means_sort.at<double>(j)) < 0.01){
                        indizes.at(i) = j;
                    }
                }
            }

            em_model->getMeans() = means_sort;

            for(size_t i = 0; i < labels.rows; i++){
                int j = labels.at<int>(i);
                labels.at<int>(i) = indizes.at(j);
            }

            cv::Mat tmp;
            em_model->getWeights().copyTo(tmp);
            
            //change order of weights
            for(size_t i = 0; i < em_model->getWeights().cols; i++){
                em_model->getWeights().at<double>(i) = tmp.at<double>(indizes.at(i));
            }
        }

        // add likelihoods to regions
        int l_i = 0; //iterator through probability msatrix
        for(size_t i = 0; i < regions.size(); i++){
            for(size_t j = 0; j < regions[i].size(); j++, l_i++){
                
                regions.at(i).at(j).likelihood.resize(n_classes,0);
                
                for(size_t k = 0; k < n_classes; k++){
                    // regions[i][j].likelihood[k] =  em_model->getPrgetProbs().at<double>(l_i, k);
                    regions[i][j].likelihood[k] = probs.at<double>(l_i, k);
                }
            
            }
        }

        Classification::printEMResults(em_model, outdir, samplevector, probs, likelihood, labels, isFlat);

    }

    /**
     * This function approximates a Gaussian Mixture Distribution (GMM) to the area of a set of regions.
     *
     *
     * The parameters of the GMM are estimated with EM using the opencv library em (http://docs.opencv.org/modules/ml/doc/expectation_maximization.html).
     * Labels are returned for classes in increasing order (first class has the smalles mean area, last class has the biggest mean area).
     * Labeling results is saved in the output directory path.
     *
     * @param K number of Gaussian functions
     * @param regions two-dimensional set of regions
     * @param labels class labels
     * @param outdir output directory path
     *
     */
    void approxMixDistrEM(int K, std::vector<std::vector<Region>> &regions, cv::Mat &labels, cv::Ptr<cv::ml::EM> &em_model, const std::string &outdir, bool isFlat){
    
        /// fill samplevector with observated areas
        std::vector<int> samplevector;
        for(size_t i = 0; i < regions.size(); i++){
            for(size_t j = 0; j < regions.at(i).size(); j++){
                int  a = regions.at(i).at(j).getArea();
                samplevector.push_back(a);
            }
        }

        int nsamples = samplevector.size();
        /// one-dimensional samples
        cv::Mat samples( nsamples, 1, CV_32FC1); 

        for(int i = 0; i < nsamples; i++){
            samples.at<float>(i) = (float) samplevector.at(i);
        }

        /// if em_model is empty, create a new object
        if ( em_model.empty() ){
            em_model = cv::ml::EM::create();
        }
        
        /// declare EM object and initialize parameter for EM object
        em_model->setClustersNumber(K); 
        em_model->setCovarianceMatrixType(cv::ml::EM::COV_MAT_SPHERICAL);
        em_model->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER+cv::TermCriteria::EPS, 500, 0.1));

        /// output of EM-model
        cv::Mat probs;
        cv::Mat likelihood;

        /// learn gmm-model
        em_model->trainEM( samples, likelihood, labels, probs );
        
        /// determine if labels are in the right order (0,1,2,...)
        cv::Mat means = em_model->getMeans();
        cv::Mat means_sort;
        cv::sort(means, means_sort, cv::SORT_EVERY_COLUMN + cv::SORT_ASCENDING);

        /// check if classes are in an ascending order
        cv::Mat means_sub;
        cv::subtract(means, means_sort, means_sub);
        means_sub = cv::abs(means_sub);

        double sum = 0;
        for(int i = 0; i < K; i++){
            sum += means_sub.at<double>(i,0);
        }

        /// change labeling if it is not in an increasing order
        if(sum > 0.1){
            
            /// lookuptable of indizes
            std::vector<int> indizes(K,0);

            /// fill lookup-table
            for(int i = 0; i < K; i++){
                double temp = means.at<double>(i);
                for(int j = 0; j < K; j++){
                    if(abs(temp - means_sort.at<double>(j)) < 0.01){
                        indizes.at(i) = j;
                    }
                }
            }
            em_model->getMeans() = means_sort;

            for(int i = 0; i < labels.rows; i++){
                int j = labels.at<int>(i);
                labels.at<int>(i) = indizes.at(j);
            }
    
            cv::Mat tmp;
            em_model->getWeights().copyTo(tmp);
                    
            /// change order of weights
            for(int i = 0; i < em_model->getWeights().cols; i++){
                em_model->getWeights().at<double>(i) = tmp.at<double>(indizes.at(i));
            }
    
        }

        /// add likelihoods to regions
        int l_i = 0; //iterator through probability msatrix
        for(size_t i = 0; i < regions.size(); i++){
            for(size_t j = 0; j < regions.at(i).size(); j++, l_i++){
                regions.at(i).at(j).likelihood.resize(K,0);
                for(int k = 0; k < K; k++){
                    regions.at(i).at(j).likelihood.at(k) = probs.at<double>(l_i, k);
                }
            }
        }

        Classification::printEMResults(em_model, outdir, samplevector, probs, likelihood, labels, isFlat);

    }

    void approxMixDistrEM(int K, std::vector<std::vector<Region*>> &regions, cv::Mat &labels, cv::Ptr<cv::ml::EM> &em_model, const std::string &outdir, bool isFlat){

        /// fill samplevector with observated areas
        std::vector<int> samplevector;
        for(size_t i = 0; i < regions.size(); i++){
            for(size_t j = 0; j < regions.at(i).size(); j++){
                int  a = regions.at(i).at(j)->getArea();
                samplevector.push_back(a);
            }
        }

        int nsamples = samplevector.size();
        /// one-dimensional samples
        cv::Mat samples( nsamples, 1, CV_32FC1); 

        for(int i = 0; i < nsamples; i++){
            samples.at<float>(i) = (float) samplevector.at(i);
        }

        /// if em_model is empty, create a new object
        if ( em_model.empty() ){
            em_model = cv::ml::EM::create();
        }

        /// declare EM object and initialize parameter for EM object
        em_model->setClustersNumber(K); 
        em_model->setCovarianceMatrixType(cv::ml::EM::COV_MAT_SPHERICAL);
        em_model->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER+cv::TermCriteria::EPS, 500, 0.1));
        
        /// output of EM-model
        cv::Mat probs;
        cv::Mat likelihood;

        /// learn gmm-model
        em_model->trainEM( samples, likelihood, labels, probs );

        /// determine if labels are in the right order (0,1,2,...)
        cv::Mat means = em_model->getMeans();
        cv::Mat means_sort;
        cv::sort(means, means_sort, cv::SORT_EVERY_COLUMN + cv::SORT_ASCENDING);
    
        /// check if classes are in an ascending order
        cv::Mat means_sub;
        cv::subtract(means, means_sort, means_sub);
        means_sub = cv::abs(means_sub);

        double sum = 0;
        for(int i = 0; i < K; i++){
            sum += means_sub.at<double>(i,0);
        }

        /// change labeling if it is not in an increasing order
        if(sum > 0.1){
            
            /// lookuptable of indizes
            std::vector<int> indizes(K,0); 

            /// fill lookup-table
            for(int i = 0; i < K; i++){
                
                double temp = means.at<double>(i);
                
                for(int j = 0; j < K; j++){
                    if(abs(temp - means_sort.at<double>(j)) < 0.01){
                        indizes.at(i) = j;
                    }
                }
            }
            em_model->getMeans() = means_sort;

            for(int i = 0; i < labels.rows; i++){
                int j = labels.at<int>(i);
                labels.at<int>(i) = indizes.at(j);
            }
    
            cv::Mat tmp;
            em_model->getWeights().copyTo(tmp);
            
            /// change order of weights
            for(int i = 0; i < em_model->getWeights().cols; i++){
                em_model->getWeights().at<double>(i) = tmp.at<double>(indizes.at(i));
            }
        
        }

        /// add likelihoods to regions
        int l_i = 0; //iterator through probability msatrix
        for(size_t i = 0; i < regions.size(); i++){
            for(size_t j = 0; j < regions.at(i).size(); j++, l_i++){
                regions.at(i).at(j)->likelihood.resize(K,0);
                for(int k = 0; k < K; k++){
                    regions.at(i).at(j)->likelihood.at(k) = probs.at<double>(l_i, k);
                }
            }
        }

        Classification::printEMResults(em_model, outdir, samplevector, probs, likelihood, labels, isFlat);

    }

    /**
     * This function prints the results of GMM classification to an output file.
     *
     * @param em_model EM-Model
     * @param outdir output file
     * @param samplevector set of samples
     * @param labels classification labels
     *
     */
    void printEMResults(cv::Ptr<cv::ml::EM> &em_model, const std::string &outdir, const std::vector<int> &samplevector, cv::Mat &probs, cv::Mat &likelihood, cv::Mat &labels, bool isFlat) {

        // print results of GMM to file
        ofstream outfile;
        std::string file = outdir + "em.txt";
        
        // io::create_directory(file);
        outfile.open(file.c_str());
        
        if(outfile.is_open()){

            outfile << "means: " << std::endl;
            outfile << em_model->getMeans() << std::endl;

            outfile << "variances: " << std::endl;
            std::vector<cv::Mat> covariances;
            em_model->getCovs(covariances);
            
            for(size_t i = 0; i < covariances.size(); i++)
                outfile << covariances.at(i) << std::endl;

            outfile << "weights: " << std::endl;
            outfile << em_model->getWeights() << std::endl;

            outfile << "probabilities: " << std::endl;
            outfile << probs << std::endl;

            outfile << "Likelihood: " << std::endl;
            outfile << likelihood << std::endl;

            outfile << "Labeling: " << std::endl;
            outfile << labels << std::endl;

            outfile.close();
        }


        // print likelihoods
        file = outdir + "likelihoods.txt";
        
        outfile.open(file.c_str());
        int cols = probs.cols;

        if(outfile.is_open()){
            // header
            outfile << "A\t" ;
            for(size_t j = 0; j < cols; j++){
                outfile << "p" << j << "\t";
            }
            outfile << "label" << std::endl;

            // data
            for(size_t i = 0 ; i< samplevector.size(); i++){
                outfile << samplevector.at(i) << "\t";

                for(size_t j = 0; j < cols; j++){
                    outfile << probs.at<double>(i,j) << "\t";
                }

                outfile << labels.at<int>(i) << std::endl;
            }
            outfile.close();
        }

    }	


} // Classification