/*  
Copyright by Jan-Philipp_Praetorius

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

#include "segment_via_gmmFCs.h"
#include "segment_via_GMM.h"
#include "IOputs.h"
#include "Classification.h"
#include "Segmentation.h"
#include "RegionSplitting.h"


namespace gmmFCs 
{

    /**
     * This function reads the binary result images of the GMM classification method
     * then saves all images and regions in the corresponding vectors
     *
     * @param outputdir directory to the binary image files
     * @param images empty vector; to be filled with images from the directory
     * @param regions empty vector of vectors of regions; to be filled with regions from the images
     */
    void fillBR(const std::string outputdir, std::vector<cv::Mat> &images, std::vector<std::vector<Region>> &regions, int &count, const bool &clear_border){
        std::vector<Region> binary;
        
        regions.clear();
        images.clear();
        count = 0;

        io::read_images(outputdir, images, cv::IMREAD_GRAYSCALE);

        for(size_t i = 0; i < images.size(); ++i){
            Segmentation::segmentRegionPixelBased2DIterative(images[i], binary, clear_border);
            regions.push_back(binary);
            
            count += binary.size();
        }
        
    }

    /**
     * Go through all regions and save regions, that do not overlap with fungal regions or dead cell regions
     */
    void classify_PMN(const int &n_classes, std::vector<std::vector<Region>> &binaryRegions, cv::Mat &labels, const std::string &OUTPUTDIR, double &d){

        std::string OUT_EM = io::create_directory(OUTPUTDIR, "em/");

        // EM-model and output (label) of the model
        cv::Ptr<cv::ml::EM> em_model = cv::ml::EM::create();
    
        Classification::approxMixDistrEM(n_classes, binaryRegions, labels, em_model, OUT_EM, true);

		int counter = 0;
		int counterArea = 0;
		std::cout << "Calculate average diameter of PMN regions..." << std::endl;

        // calculate average diameter of the cells; exclude Objects of class 0 (= Noise)
        
        ///// notes from susanne: average diameter will be used later?!  /////
		// is used in the combineRegionsOverlap function (which is called by addingFlat)
		
		for(size_t i = 0; i < binaryRegions.size(); i++){ 
			for(size_t j = 0; j < binaryRegions.at(i).size(); j++){

                cell_class::immune mylabel = cell_class::setCellType(labels.at<int>(counter));
				
                binaryRegions.at(i).at(j).setClass( mylabel );
				
                if(binaryRegions.at(i).at(j).getClass() == cell_class::immune::SINGLE ){
					d += 2*sqrt(binaryRegions.at(i).at(j).getArea()/M_PI);
					counterArea++;
				}
				
                counter++;
			}
		}

        d = d / counterArea;

		std::cout << "Avg. diameter " << d << std::endl;
		std::cout << "Counter area " << counterArea << std::endl;
		std::cout << "Counter " << counter << std::endl;
        std::cout << "Classify all regions according to size..." << std::endl;
		
        labels.resize(0);
		
        OUT_EM = io::create_directory(OUTPUTDIR, "em_sizes/");
        Classification::approxMixDistrEM(n_classes, binaryRegions, labels, em_model, OUT_EM, true);
		
        counter = 0;
		// set class of regions
		for(int unsigned i = 0; i < binaryRegions.size(); i++){
			for(int unsigned j = 0; j < binaryRegions.at(i).size(); j++){

                cell_class::immune mylabel = cell_class::setCellType(labels.at<int>(counter));

				binaryRegions.at(i).at(j).setClass( mylabel );
				counter++;
			}
		}

    }

    /**
     *
     * @param regions binary Regions of cells
     * @param rows rows of output image
     * @param cols columns of output image
     * @param labels labels of the regions infered through gmm
     * @param outputdir directory to save the generated images showing the classified cells in different grey values
     *
     */
    void create_gmm_image(std::vector<std::vector<Region>> &regions, const int &rows, const int &cols, const cv::Mat &labels, const std::string &outputdir, const int &K, const int &N_TEMPVAR) {
        
        std::vector<cv::Vec2i>::iterator vecit;
        std::vector<std::vector <Region> >::iterator vreg_it;
        std::vector<Region>::iterator reg_it;
    
        int i = 0, j = 0;
        for (vreg_it = regions.begin(); vreg_it != regions.end(); vreg_it++){
            
            cv::Mat tmp(rows, cols, CV_8UC1, cv::Scalar(0));
            cv::Mat im2(rows, cols, CV_8UC1, cv::Scalar(0));
            
            for (reg_it = vreg_it->begin(); reg_it != vreg_it->end(); reg_it++){

                int l = labels.at<int>(i);

                for(vecit = reg_it->region_pixels.begin(); vecit != reg_it->region_pixels.end(); vecit++){
                    cv::Point p(vecit->val[1], vecit->val[0]);

                    for (int k = 0; k < K; k++){
                        if (l == k){
                            tmp.at<uchar>(p) = floor((k+1)*(255/K));
                        }
                    }
                }
                i++;
            }
        
            cv::imwrite(outputdir + std::to_string(j+ int(floor(N_TEMPVAR/2)))+"classification_by_size.png", tmp);
            j++;

	    }

    }

    /**
     * This function applies morphological operations (opening then closing)
     * using a kernel as a limit.
     *
     * @param m an image contains contains the result from the combination of
     * static and mobile elements(GMM)
     * @param kernel the kernel used in the morph. ops
     */
    void enhncedOpen (cv::Mat &m, const cv::Mat &kernel){
        // Opening
        cv::Mat temp;
        m.copyTo(temp);
        cv::erode(temp,m,kernel,cv::Point(-1,-1),1,cv::BORDER_CONSTANT,cv::morphologyDefaultBorderValue());
        m.copyTo(temp);
        cv::dilate(temp,m,kernel,cv::Point(-1,-1),1,cv::BORDER_CONSTANT,cv::morphologyDefaultBorderValue());
    }

    void enhncedClose (cv::Mat &m, const cv::Mat &kernel){
        // Closing
        cv::Mat temp;
        m.copyTo(temp);
        cv::dilate(temp,m,kernel,cv::Point(-1,-1),1,cv::BORDER_CONSTANT,cv::morphologyDefaultBorderValue());
        m.copyTo(temp);
        cv::erode(temp,m,kernel,cv::Point(-1,-1),1,cv::BORDER_CONSTANT,cv::morphologyDefaultBorderValue());
    }

    /**
     * This function sets the region's class (noise, single cell, cell cluster)
     * into binaryRegions variable(vector of regions) from labels (Single array matrix)
     *
     * @param regions vector of regions needs to be checked for noise
     * @param rows number of rows in the image set
     * @param cols number of columns in the image set
     *
     *
     * - collect areas of the regions
     * - search for the class with the smallest mean value in the em model
     * - go through all regions and erase regions which are predicted to have that class (small area --> noise class)
     * - I think this will lead to the erasure of mostly fungal cells. which doesn't matter, because this function is called to erase noise in the images
     *   that are used to add the flat cells. --> Therefore this procedure will also prevent the adding of too small cells
     */
    void removeNoiseRegions(cv::Ptr<cv::ml::EM> em_model, std::vector <Region> &regions, const int &rows, const int &cols) {

        std::vector<cv::Point> points;
        cv::RotatedRect rect;
        float dmax;
        std::vector<int> sizeVect;
        
        for(size_t i = 0; i < regions.size(); i++){
            sizeVect.push_back(regions.at(i).getArea());
        }

        // int n = sizeVect.size();        
        // cv::Mat samples( 1, 1, CV_32FC1); //one-dimensional samples
        // cv::Mat out( n, 1, CV_32FC1); //one-dimensional samples
        
        cv::Vec2d response; // = cvRound ( em_model.predict ( samples., NULL ) );
        int noiseClass = 0; //= em_model.get_means()->rows;
        
        cv::Mat means = em_model->getMeans();
        std::vector<int> meansList = means.reshape(0, 1);
        int temp = meansList.at(0);
        
        for(size_t i = 0; i< meansList.size() ;i++){
            if(meansList.at(i) < temp){
                noiseClass = i;
                temp = meansList.at(i);
            }
        }

        // CAUTION: in samples just one feature will be added, but model needs 2-dim features
        cv::Mat samples = cv::Mat( 1 , means.cols, CV_32FC1); 
        // cv::Mat out( n, 1, CV_32FC1); //one-dimensional samples

        for(size_t i = 0; i < regions.size(); i++){
            samples.at<float>(0) = (float) sizeVect.at(i);
            response =  em_model->predict ( samples );
            
            if(response.val[0] == noiseClass){
                regions.erase(regions.begin()+i);
                sizeVect.erase(sizeVect.begin()+i); 
                
                --i;
            }
        }
 
        for(size_t i = 0; i < regions.size(); i++){
            for(size_t ci = 0; ci < regions.at(i).region_pixels.size(); ci++){
                cv::Point p(regions.at(i).region_pixels.at(ci)[1], regions.at(i).region_pixels.at(ci)[0]);
                points.push_back(p);
            }

            // try-catch-clause added by Philipp Praetorius to prevent
            //  OpenCV Error: Incorrect size of input array (There should be at least 5 points to fit the ellipse) in fitEllipse
            try
            {
                rect = cv::fitEllipse(points);
            }
            catch(const std::exception& e)
            {
                std::cout << "cv::fitEllipse: There should be at least 5 points to fit the ellipse in function fitEllipse: " << points.size() << " Points" << std::endl;
                std::cerr << e.what() << '\n';

                points.clear();
                std::cout << "CAUTION: clear points-vector and do not erase any region in regions-vector" << std::endl;
                continue;                
            }

            if(rect.size.height > rect.size.width )
                dmax = rect.size.height;            
            else
                dmax = rect.size.width;
            
            if((int)dmax >= (rows+cols)/8 || gmm::is_line(regions.at(i))){
                regions.erase(regions.begin()+i);
                --i;
            }
            
            points.clear();
        }
            
    }

    /**
     * written by Stefanie
     *
     * @param output vector of regions return the final result
     * @param binary vector of regions contains the existing result regions
     * @param classified vector of regions contains the result from the combination of static and mobile elements(GMM) after removing noise
     * @param added vector of regions contains all newly added regions
     *
     * 1) go through all FC candidate regions
     * 		- go through all already found regions
     * 			- if the candidate and one of the already found regions overlap:
     * 				- check, if the distance of the centroids of the regions is less than d
     * 					- calculate the amount of overlap --> if the old region completely resides in the new region:
     * 						- check, if the FC candidate is 4 times larger than the already found region
     * 							- add the candidate
     * 						- try to split the region into 2 subregions using ellipse fitting
     * 						- erase the originl region (?)
     * 						- go through the subregions:
     * 							- compute overlap of the subregion and the already found region
     * 								- if the already found region is engulfed --> remove the new subregion
     * 						- go through all remaining subregions:
     * 							- save them in the FC-candidate vector
     * 						- start from the top of the candidate vector
     * 			- if no overlap was found --> save the candidate
     *		- if there are no previous found regions --> save the candidate
    *
    *
    */
    void combineRegionsOverlap(std::vector<Region> &output, std::vector<Region> &allRegions, std::vector<Region> &classified, std::vector<Region> &morph, std::vector<Region> &added, const double &d) {
        
        for(size_t i = 0; i < allRegions.size();i++){
            output.push_back(allRegions.at(i));
        }

        for(size_t i = 0; i < classified.size();i++){
            if(allRegions.size() != 0){
                for (size_t j = 0; j < allRegions.size();j++){
                    
                    if(overlap(classified.at(i), allRegions.at(j))){
                        if(cv::norm(classified.at(i).getCentroidPoint() - allRegions.at(j).getCentroidPoint()) < d){
                            
                            if(overlapN(classified.at(i),allRegions.at(j)) == allRegions.at(j).getArea()){

                                // existing area is smaller than a fourth of the new area --> add new area directly
                                if(allRegions.at(j).getArea() < 0.25*classified.at(i).getArea()){ 
                                    added.push_back(classified.at(i));
                                }

                                std::vector<Region> subregions(2);

                                RegionSplitting::separateClusterGMM(classified.at(i), subregions);

                                classified.erase(classified.begin()+i);
                                if (subregions.size() != 0){
                                    for(size_t k = 0; k < subregions.size();k++){
                                        if(overlapN(subregions.at(k),allRegions.at(j)) == allRegions.at(j).getArea()){
                                            subregions.erase(subregions.begin()+k);
                                        }
                                    }
                                }

                                if (subregions.size() != 0){
                                    for(size_t k = 0; k < subregions.size();k++){
                                        classified.push_back(subregions.at(k));
                                    }
                                }
                                i = 0;
                            }
                        }
                        break;
                    }

                    if(j == allRegions.size()-1){
                        output.push_back(classified.at(i));
                        added.push_back(classified.at(i));
                    }
                }
            }
            else{
                output.push_back(classified.at(i));
                added.push_back(classified.at(i));
            }
        }

        for(size_t i = 0; i < morph.size();i++)
            for (size_t j = 0; j < output.size();j++){
                if(overlap(morph.at(i),output.at(j))){
                    break;
                }
                if(j == output.size()-1){
                    output.push_back(morph.at(i));
                    added.push_back(morph.at(i));
                }
            }

    }

    /**
    * This function draws an image of binary regions.
    *
    * @param im output image
    * @param regions vector of regions needs to be drawn
    */
    void drawImage(cv::Mat & im, std::vector<Region> & regions) {
        std::vector<cv::Vec2i>::iterator vecit;
        for(size_t i = 0; i<regions.size();i++){
            for(vecit =  regions.at(i).region_pixels.begin(); vecit != regions.at(i).region_pixels.end(); vecit++){
                cv::Point p(vecit->val[1], vecit->val[0]);
                im.at<uchar>(p) = 255;
            }
        }
    }    

    /**
     *
     * @param c image containing result from the GMM classification
     * @param binary already segmented regions
     * @param fungalbinary fungal regions
     * @param deadcellbinary dead cell regions
     * @param add empty image for saving flat cells
     * @param imgnum number of current image
     *
     * @return result image with all cells
     *
     * 1. Threshold the GMM image at 126 and clone the result
     *    - save all found regions in the classf region vector
     * 2. enhance the cloned image with the prepared kernel (perform one opening and one closing operation)
     *    - save all resulting regions in the morph region vector
     *
     */
    cv::Mat addingFlat(cv::Ptr<cv::ml::EM> em_model, cv::Mat c, std::vector<Region> binary, std::vector<Region> fungalbinary, std::vector<Region> deadcellbinary, cv::Mat &add, const double &d, const bool &clear_border){

        std::vector<Region> morph;
        std::vector<Region> classf;
        std::vector<Region> output;
        std::vector<Region> added;

        // generate empty images
        cv::Mat m;
        cv::Mat temp;
        cv::Mat result = cv::Mat::zeros(c.rows,c.cols,0);
        add = cv::Mat::zeros(c.rows,c.cols,0);

        // generate kernel for morphological operation
        cv::Mat kernel = cv::getStructuringElement(cv::MORPH_OPEN,cv::Size(15,15),cv::Point(7,7)); //default size (15,15)

        // c is the GMM image after thresholding by 126
        c.copyTo(temp);
        cv::threshold(temp,c,126,255,cv::THRESH_BINARY);
        m = c.clone();
        
        Segmentation::segmentRegionPixelBased2DIterative(c, classf, clear_border);

        // m is the thresholded image after a morphological operation
        gmmFCs::enhncedOpen(m, kernel);
        gmmFCs::enhncedClose(m, kernel);

        Segmentation::segmentRegionPixelBased2DIterative(m, morph, clear_border);

        cv::Mat tmp = gmm::create_binary_image_from_regions(c.rows, c.cols, binary);        

        // erase regions of class 0 (=Noise) from already found regions only erase them, 
        // if they don't overlap with fungal regions
        // use classification already done in the main function before calling GMM_from_thread_FC
        for(size_t i = 0; i < binary.size(); i++){
            
            if(binary.at(i).getClass() == cell_class::immune::NOISE ){
                
                std::vector<int> ovlp;
                for (size_t j = 0; j < fungalbinary.size(); j++){
                    int o = overlapN(binary.at(i), fungalbinary.at(j));
                    if (o > 0){
                        ovlp.push_back(o);
                    }
                }
                
                for (size_t j = 0; j < deadcellbinary.size(); j++){
                    int o = overlapN(binary.at(i), deadcellbinary.at(j));
                    if (o > 0){
                        ovlp.push_back(o);
                    }
                }
                
                if (ovlp.size() == 0){
                    binary.erase(binary.begin()+i);
                    --i;
                }
            }
       
        }
        
        tmp = gmm::create_binary_image_from_regions(c.rows, c.cols, binary);
        
        // after removing noise, add all the fungal and dead cell regions
        cv::Mat tmp2 = gmm::create_binary_image_from_regions(c.rows, c.cols, fungalbinary);
        cv::Mat tmp3 = gmm::create_binary_image_from_regions(c.rows, c.cols, deadcellbinary);
        cv::add(tmp, tmp2, tmp);
        cv::add(tmp, tmp3, tmp);

        std::vector<Region> allRegions;
        Segmentation::segmentRegionPixelBased2DIterative(tmp, allRegions, clear_border);

        // erase regions with area <= 5
        int count = 0;
        for(size_t i = 0; i < classf.size(); i++){
            if (classf.at(i).getArea() <= 5){
                classf.erase(classf.begin()+i);
                count ++;
                --i;
            }
        }
        
        tmp = gmm::create_binary_image_from_regions(c.rows, c.cols, classf);
        
        // erase regions with area <= 5
        count = 0;
        for(size_t i = 0; i < morph.size(); i++){
            if (morph.at(i).getArea() <= 5){
                morph.erase(morph.begin()+i);
                count ++;
                --i;
            }
        }
        
        tmp = gmm::create_binary_image_from_regions(c.rows, c.cols, morph);
        
        // erase regions that are very elongated (lines)
        count = 0;
        for(size_t i = 0; i < classf.size(); i++){
            if(gmm::is_line(classf.at(i))){
                classf.erase(classf.begin()+i);
                count ++;
                --i;
            }
        }
        
        tmp = gmm::create_binary_image_from_regions(c.rows, c.cols, classf);
        
        gmmFCs::removeNoiseRegions(em_model, classf, c.rows, c.cols);
        
        tmp = gmm::create_binary_image_from_regions(c.rows, c.cols, classf);
        
        gmmFCs::removeNoiseRegions(em_model, morph, c.rows, c.cols);
        
        tmp = gmm::create_binary_image_from_regions(c.rows, c.cols, morph);
        
        // filling holes
        Segmentation::fill_holes(classf, clear_border);
        tmp = gmm::create_binary_image_from_regions(c.rows, c.cols, classf);
        
        kernel = cv::getStructuringElement(cv::MORPH_OPEN, cv::Size(5,5), cv::Point(2,2));
        gmmFCs::enhncedOpen(tmp, kernel);
        gmmFCs::enhncedClose(tmp, kernel);
        
        Segmentation::segmentRegionPixelBased2DIterative(tmp, classf, clear_border);
 
        // remove noise in classf and again remove elongated objects
        gmmFCs::removeNoiseRegions(em_model, classf, c.rows, c.cols);
        
        count = 0;
        for(size_t i = 0; i < classf.size(); i++){
            
            if(gmm::is_line(classf.at(i))){
                classf.erase(classf.begin()+i);
                count ++;
                --i;
            }
        }
 
        tmp = gmm::create_binary_image_from_regions(c.rows, c.cols, classf);
    
        Segmentation::segmentRegionPixelBased2DIterative(tmp, classf, clear_border);

        // remove regions that touch the image border
        int i = 0;
        while(i < int(classf.size())){
            if(Segmentation::touchesBorder(classf.at(i).contour_pixels, tmp)){
                classf.erase(classf.begin()+i);
                i--;
            }
            i++;
        }

        tmp = gmm::create_binary_image_from_regions(c.rows, c.cols, classf);
        
        // compensate oversegmentation by erosion
        kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(9,9), cv::Point(4,4));
        tmp2 = gmm::create_binary_image_from_regions(c.rows, c.cols, classf);
        cv::erode(tmp2, tmp, kernel, cv::Point(-1,-1), 1, cv::BORDER_CONSTANT, cv::morphologyDefaultBorderValue());
        
        Segmentation::segmentRegionPixelBased2DIterative(tmp, classf, clear_border);
 
        // start combining Regions
        gmmFCs::combineRegionsOverlap(output, allRegions, classf, morph, added, d);
 
        for(size_t i = 0; i < added.size(); i++){
            if(added.at(i).getArea() < 50){
                added.erase(added.begin()+i);
                --i;
            }
        }
        
        // for all
        drawImage(result,output);
        // for flat cells only
        drawImage(add,added);
        
        // this opening operation was added
        kernel = cv::getStructuringElement(cv::MORPH_OPEN, cv::Size(11,11), cv::Point(5,5));
        gmmFCs::enhncedOpen(add, kernel);

        return result;
    }

    /**
     * General function for GMM-segmentation adding flat cells.
     *
     * @param tid 		thread ID
     * @param inputdir 	input directory path
     * @param outputdir output directory path
     * @param outputdirGMM output directory path for gmm images
     */
    void perform_gmm_segmentation(const std::vector<cv::Mat> &images, const int &N_THREADS, const int &N_TEMPVAR, const int &C, const double &d,
            std::vector<std::vector<Region>> binaryRegions, std::vector<std::vector<Region>> fungalbinaryRegions, std::vector<std::vector<Region>> deadcellbinaryRegions, 
            const std::string OUTPUTDIR, const bool &clear_border, const bool &FLAG_DEBUG) {

        std::string OUTPUTDIR_GMM = io::create_directory( OUTPUTDIR, "gmm/" );
        std::string OUTPUTDIR_FC = io::create_directory( OUTPUTDIR, "binaryFC/" );

        int first = int(floor(N_TEMPVAR/2));
        int last = (signed) images.size() - first -1;

        // for all images allocated to this thread (excluding the very first and very last images because of the delay), do:
        // multi-threaded processing of brightfield images
        #pragma omp parallel for firstprivate(images) num_threads(N_THREADS)
        for(size_t i = 0; i <= images.size(); ++i) {

            // compute gmm of the second to penultimate image
            if(i >= (size_t)first && i <= (size_t)last){
                cv::Mat img_gmm; /*! image with pixels colored in different grey values according to the corresponding gmm class */
                cv::Mat img_binary; /*! new binary images containing all cells */
                cv::Mat img_added;/*! binary image containing flat cells */
                
                // classify image pixels into three classes (according to spatial and temporal variance)
                cv::Ptr<cv::ml::EM> em_model;
                gmm::classify_pixel(i, img_gmm, images, em_model, N_TEMPVAR, C, FLAG_DEBUG);
                
                // save FC image
                io::write_image(OUTPUTDIR_GMM, img_gmm, i, true);
               
                std::cout << "Number of Binary Objects = "<< binaryRegions.at(i-first).size() << std::endl;
                
                img_binary = addingFlat(em_model, img_gmm, binaryRegions.at(i-first), fungalbinaryRegions.at(i), deadcellbinaryRegions.at(i), img_added, d, clear_border);
                            
                std::vector<std::string> names = {"00_c","01_m","02_bin","03_bin","04_c","05_m","06_c","07_c","08_m","09_c", "10_c", "11_c", "12_c", "13_c", "14_"};
                
                // save FC image
                io::write_image(OUTPUTDIR_FC, img_added, i, true);
                 
            }
        
        }

    }

    /**
     * Main-Function for GMM-segmentation
     *
     * @param ...
     */
    void gmmFCs_segmentation(const int &n_classes, const std::vector<cv::Mat> &images, std::vector<cv::Mat> &images_fungi, const int &N_THREADS, const int &N_TEMPVAR, const std::string &INPUTDIR, const std::string &OUTPUTDIR, const bool &clear_border, const bool &FLAG_DEBUG){

        std::vector<cv::Mat> images_binary, images_dead;

        ///// fill region vector
        int n_binary, n_fungi, n_dead;
        std::vector<std::vector<Region>> binaryRegions, fungalbinaryRegions, deadcellbinaryRegions;
        // CAUTION: self-made in/output-path
        std::string IN_BINARY = INPUTDIR + "binary/";
        std::string IN_GREEN = INPUTDIR + "green/";
        std::string IN_RED = INPUTDIR + "red/";
        
        std::string OUT_CLASS = OUTPUTDIR + "classificationBySize/";
        

        gmmFCs::fillBR(IN_BINARY, images_binary, binaryRegions, n_binary, clear_border);
        std::cout << " [0.5] Loading segmented fungal regions ..." << std::endl;
		gmmFCs::fillBR(IN_GREEN, images_fungi, fungalbinaryRegions, n_fungi, clear_border);
		std::cout << " [0.6] Loading segmented dead cell regions ..." << std::endl;
		gmmFCs::fillBR(IN_RED, images_dead, deadcellbinaryRegions, n_dead, clear_border);

        std::cout << "Binary images: " << binaryRegions.size() << std::endl;
		std::cout << "Fungal images: " << fungalbinaryRegions.size() << std::endl;
		std::cout << "Dead cell images: " << deadcellbinaryRegions.size() << std::endl;

		std::cout << "Binary regions: " << n_binary << std::endl;
		std::cout << "Fungal regions: " << n_fungi << std::endl;
		std::cout << "Dead cell regions: " << n_dead << std::endl;

		std::cout << "Classify only PMN regions according to size..." << std::endl;
		std::cout << "\t - collect PMN regions" << std::endl;

        std::cout << "\t - classify binary regions according to size" << std::endl;
        
        cv::Mat labels;
        double d;
        gmmFCs::classify_PMN(n_classes, binaryRegions, labels, OUTPUTDIR, d);

		std::cout << "EM done. Output images are generated." << std::endl;
        
        gmmFCs::create_gmm_image(binaryRegions, images_binary[0].rows, images_binary[0].cols, labels, OUT_CLASS, n_classes, N_TEMPVAR);

		std::cout << "Generation of output images done." << std::endl;

		std::cout << "------------------ Start FC segmentation ------------------------" << std::endl;

        gmmFCs::perform_gmm_segmentation(images, N_THREADS, N_TEMPVAR, n_classes, d, binaryRegions, fungalbinaryRegions, deadcellbinaryRegions, OUTPUTDIR, clear_border, FLAG_DEBUG);

		std::cout<<"Segmentation of flat cells done." << std::endl;

    }

} // gmmFCs