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

#include <iostream>
#include <sstream>
#include <chrono>
#include <fstream>
#include <dirent.h>
#include <errno.h>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/program_options.hpp>

#include "IOputs.h"
#include "ImageProcessingToolbox.h"

namespace IPT = ImageProcessingToolbox;


namespace io {

	/**
	 * Search for specified file in specified directory
	 * 
	 * @param dirName string-path of the directory 
	 *
	 */
	boost::optional<path> find_file(const path& dir_path, const path& file_name) {
  		const recursive_directory_iterator end;
  		const auto it = std::find_if(recursive_directory_iterator(dir_path), end,
                          [&file_name](const directory_entry& e) {
                            return e.path().filename() == file_name;
                          });
  		return it == end ? boost::optional<path>() : it->path();
	}
		
	/**
	 * Create specified directory with optional suffix parameter with respect to the last char
	 * 
	 * @param dirName string-path of the directory 
	 * 
	 * @return path of created directory
	 *
	 */
	std::string create_directory(const std::string &dirName, const std::optional<std::string> suffix){

		std::stringstream stringStream;
		stringStream << dirName;

		if( suffix.has_value() )  {
		
			char last_char = dirName.back();

			// check the last char for a slash to ensure store it in the outputfolder
			if (last_char != '/') {
				stringStream << "/";
			}

			stringStream << suffix.value();
		}
		
		boost::filesystem::path dir(stringStream.str());

		if (boost::filesystem::exists(dir)) { 
			std::cout << " [INFO] IOputs::createDirectory: " << stringStream.str() << " already exist" << std::endl; 
		} else {
			bool res = boost::filesystem::create_directories(dir);
			std::cout << " [INFO] Create directory: " << stringStream.str() << "\tworked: " << res << std::endl;	
		}

		return stringStream.str();	

	}

	/**
	 * Read images from specified directory to cv::Mat variable (gray-scaled)
	 * 
	 * @param dir corresponding directory
	 * @param images pointer on images-variable 
	 *
	 */
	void read_images(const std::string &dir, std::vector<cv::Mat> &images, cv::ImreadModes COLOR_FLAG, std::optional<int> num_images ) {
		
		std::cerr << " [INFO] OpenCV IMREAD-mode: " << IPT::type2str(COLOR_FLAG) << "\tint: " << COLOR_FLAG << std::endl;		

		cv::Size mysize = cv::Size(-1,-1);

		std::vector<std::string> files;
		
		// extract filename of specified directory
		io::read_filenames(dir, files);

		// for specified n_I: read specified number of images
		for(size_t i = 0; i < files.size() && i < (size_t)num_images.value_or(files.size()); ++i) { 
			
			cv::Mat img, img_cvt;

			try{  

				/// convert to 8UC1 image for imread-mode: IMREAD_GRAYSCALE
				if (COLOR_FLAG == cv::IMREAD_GRAYSCALE){

					img = cv::imread(files[i], cv::IMREAD_ANYDEPTH);
					
					double  minVal,  maxVal;
					cv::minMaxLoc(img, &minVal, &maxVal); 	

					img.convertTo(img_cvt, CV_8UC1, (1.0 / maxVal)*255. , 0. );	

				} 
				/// convert to 8UC3 image for imread-mode: IMREAD_COLOR
				else if (COLOR_FLAG == cv::IMREAD_COLOR){

					img = cv::imread(files[i], cv::IMREAD_ANYCOLOR);
					
					double  minVal,  maxVal;
					cv::minMaxLoc(img, &minVal, &maxVal); 	

					img.convertTo(img_cvt, CV_8UC3, (1.0 / maxVal)*255. , 0. );	

				}				
				else {
					img = cv::imread(files[i], cv::IMREAD_ANYCOLOR);
					img.copyTo(img_cvt);
				}	

				/// check for an empty image, skip (alternative: add empty(zero) ) image
				if ( img_cvt.empty() ) {
					std::cerr << " [WARNING] io::read_images - empty image with size: " << img_cvt.size() << "\tindex: " << i << std::endl;
					
					/// add empty image with respect to the specified COLOR-mode
					if (COLOR_FLAG == cv::IMREAD_GRAYSCALE) {
						img_cvt = cv::Mat::zeros(mysize, CV_8UC1);
					} else if (COLOR_FLAG == cv::IMREAD_COLOR) {
						img_cvt = cv::Mat::zeros(mysize, CV_8UC3);
					} else {
						img_cvt = cv::Mat::zeros(mysize, cv::IMREAD_ANYCOLOR);
					}
					
					/// skip image (might causes inconsistency when more channels will be segmented)
					// continue;
				}

				/// extract just on time the image size
				if(mysize.height == -1 and mysize.width == -1){
					mysize = img_cvt.size();
					std::cerr << " [INFO] io::read_images - image size: " << mysize << std::endl;		
				}

				images.push_back(img_cvt);
				
			}
			catch(const std::exception& e) {
				std::cerr << e.what() << std::endl;
				std::cerr << "io::read_images: " << files[i] << std::endl;										
			}

		}
	}

	/**
	 * Just read the path of the images from specified directory and store them in files
	 * 
	 * @param dir corresponding directory
	 * @param files pointer on images-file vector 
	 *
	 */
	void read_filenames(const std::string &dir, std::vector<std::string> &files) {

		// extract filename of specified directory
		path p(dir);
		if(is_directory(p)) {
			for(const auto& entry : boost::make_iterator_range(directory_iterator(p), {}))			
				files.push_back(entry.path().string());
		}
		
		// take care, that images are sorted!
		sort(files.begin(), files.end());	

	}
	

	/**
	 * Saves an image to the passed directory.
	 *
	 * @param inputdir 		is used to look up the file name of the related input image
	 * @param outputdir 	output directory
	 * @param image 	 	image to be saved
	 * @param i			 	number of the image in the file list
	 * @param addon 		(optional) suffix in filename 
	 */
	void write_image(const std::string &outDir, const cv::Mat &image, const int &i, std::optional<bool> verbose, std::optional<std::string> addon){ 
		
		std::string img_type = ".png";
		char last_char = outDir.back();

		// build string with 0-prefix for index (001, 054, ...)
		std::stringstream stringStream;

		stringStream << outDir;

		// check the last char for a slash to ensure store it in the outputfolder
		if (last_char == '/') {
			stringStream << std::setfill('0') << std::setw(3) << i;
		} else {
			stringStream << "/" << std::setfill('0') << std::setw(3) << i;
		}

		if( addon.has_value() )  {
			stringStream << "_" << addon.value();
		} 

		stringStream << img_type;

		cv::imwrite(stringStream.str(), image);

		// print save path if verbose is true
		if( verbose.has_value() and verbose.value() == true  ){
			std::cout << "save image to:\t" << stringStream.str() << std::endl;
		}		

	}

	/**
	 * Draw temp-variance, spatial-variance and histogram image for FLAG_DEBUG
	 */
	void draw_images(const cv::Mat &imtempvar, const cv::Mat &imvar){
		std::string fname = "_gmm_1_tempvar.png";

		cv::imshow("tempvar", imtempvar);
		cv::imwrite(fname, imtempvar);
		////cv::waitKey(0);
		cv::destroyAllWindows();
		cv::imshow("spatvar", imvar);
		fname = "_gmm_2_spatvar.png";
		cv::imwrite(fname, imvar);
		////cv::waitKey(0);
		cv::destroyAllWindows();

		int histSize = 256;
		float range[] = {0,256};
		const float* histRange = {range};
		bool uniform = true; bool accumulate = false;
		cv::Mat t_hist, s_hist;
		cv::calcHist(&imtempvar, 1, 0, cv::Mat(), t_hist, 1, &histSize, &histRange, uniform, accumulate);
		cv::calcHist(&imvar, 1, 0, cv::Mat(), s_hist, 1, &histSize, &histRange, uniform, accumulate);

		int hist_w = 5120; int hist_h = 4000;
		int bin_w = cvRound((double) hist_w/histSize);
		cv::Mat histImage( hist_h, hist_w, CV_8UC3, cv::Scalar( 0,0,0) );
		cv::normalize(t_hist, t_hist, 0, histImage.rows, cv::NORM_MINMAX, -1, cv::Mat() );
		cv::normalize(s_hist, s_hist, 0, histImage.rows, cv::NORM_MINMAX, -1, cv::Mat() );
		for(size_t i = 1; i < (size_t)histSize; ++i)
		{
			line( histImage, cv::Point( bin_w*(i-1), hist_h - cvRound(t_hist.at<float>(i-1)) ) ,
				cv::Point( bin_w*(i), hist_h - cvRound(t_hist.at<float>(i)) ),
				cv::Scalar( 255, 0, 0), 2, 8, 0  );
			line( histImage, cv::Point( bin_w*(i-1), hist_h - cvRound(s_hist.at<float>(i-1)) ) ,
				cv::Point( bin_w*(i), hist_h - cvRound(s_hist.at<float>(i)) ),
				cv::Scalar( 0, 255, 0), 2, 8, 0  );
		}

		/// Display
		cv::imshow("calcHist Demo", histImage );

		fname = "_hist.png";
		cv::imwrite(fname, histImage);
		cv::destroyAllWindows();
		std::cout << "\t\tsampling" << std::endl;    
	}

	/**
	 * get json file from specified path
	 */
	json get_JSON(int argc, char* argv[]){

		json J;

		std::string config_path;
    
		try {   
			namespace po = boost::program_options; 
			po::options_description desc{"Options"};
			desc.add_options()
				("config,c", po::value<std::string>(&config_path), "Directory where the config.json is located");

			po::variables_map vm;
			
			try {
				po::store(po::parse_command_line(argc, argv, desc), vm);

				if ( vm.count("help")  )
					std::cout << "Basic Command Line Parameter App: -c <path to config.json>" << std::endl << desc << std::endl; 
				po::notify(vm);
				
				if ( vm.count("config") ){
					std::cout << "Path to config.json: " << config_path << std::endl;
				} 
				else
				{
					std::cout << "Please choose path to config.json" << std::endl;
					return J;
				}            

				/// read the config.json
				try {
					/// get current path and build search-path to read config.json
					boost::filesystem::path current_path(boost::filesystem::current_path());
					const path myPath = current_path.string(); 
					const path myFile = "config.json";
					
					/// version 1: find config file in specified path
					/// config_path = io::find_file(myPath, myFile).get_value_or("not found").string();
					
					/// version 2: specify path in command line (default)
					std::ifstream i(config_path);
					
					i >> J;        
				}
				catch(const std::exception& e) {
					std::cout << " ERROR:\t JSON not found" << std::endl;
					std::cerr << e.what() << '\n';
					return 1;
				}
				
			}
			catch(po::error& e) {
				std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
				std::cerr << desc << std::endl; 
			}
	
		}
		catch (std::exception& e) {
			std::cerr << "Unhandled Exception reached the top of main: " << e.what() << ", application will now exit" << std::endl; 
		}

		return J;

	}


}
