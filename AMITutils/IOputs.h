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

#pragma once

#include <string>
#include <vector>
#include <opencv2/opencv.hpp>
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <nlohmann/json.hpp>

using namespace boost::filesystem;
using json = nlohmann::json; 


namespace io {

    int get_image_bit_depth(const int &depth);
    boost::optional<path> find_file(const path& dir_path, const path& file_name);
    std::string create_directory(const std::string &dirName, const std::optional<std::string> suffix = std::nullopt);
    void read_images(const std::string &dir, std::vector<cv::Mat> &images, cv::ImreadModes COLOR_FLAG, std::optional<int> num_images = std::nullopt);
    void read_filenames(const std::string &dir, std::vector<std::string> &files);
    void write_image(const std::string &outDir, const cv::Mat &image, const int &i, std::optional<bool> verbose = std::nullopt, std::optional<std::string> addon = std::nullopt, std::optional<std::string> image_type = std::nullopt);
    void draw_images(const cv::Mat &imtempvar, const cv::Mat &imvar);
    json get_JSON(int argc, char* argv[]);

} // io