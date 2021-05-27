

This is the first part of the algorithm for migration and interaction tracking (AMIT). In this part both immune cells and pathogens can be segmented. Multiple segmentation pipelines can be applied on multi-channel live cell videos of migration and confrontation assays (see section **Parameter instructions**).

# How to run AMIT Segmentation:

Run the executable file in the build directory by specifying the path where your configuration file with the required parameters is located with the following command: 

```console
./build/AMITSegmentation/AMITSegmentation --config <path to config file/config.json>

# example
./build/AMITSegmentation/AMITSegmentation --config ./AMITSegmentation/config.json
```

In the following table you can the see the hyper parameters for the algorithm with their corresponding meaning. You can change the values for each parameter in the **AMITSegmentation/config.json** file. To get a better understanding of the individual parameters, take a look at the section **Parameter instructions** and **Further remarks**.

|     Parameter     |      Type       |        Default         | Description                                                  |
| :---------------: | :-------------: | :--------------------: | :----------------------------------------------------------- |
|      `input`      | Directory path  |     /brightfield/      | Input directory with gray scaled images                      |
|  `input_gmmFCs`   | Directory path  |        /green/         | Input directory with green fluorescence channel  (only for `method` *gmmFCs*) |
|     `output`      | Directory path  |   /AMITSegmentation/   | Output path where all segmented images will be stored        |
|     `method`      |      Flag       |          faf           | Method used for segmentation (see **Further remarks**)       |
| `min_region_size` |   Integer > 0   |          300           | Discard objects smaller than the specified pixel size        |
| `number_closing`  |   Integer > 0   |           3            | Number of closings applied on ROI to get solid one (only for `method` *gmm*) |
| `number_tempvar`  |   Integer > 0   |           3            | Number of frames used to calculate the temporal variance (only for `method` *gmm* and *gmmFCs*) |
|   `morph_open`    |   Integer > 2   |           7            | Kernel size for morphological opening after median blur (only for `method` *canny*) |
|   `morph_close`   |   Integer > 2   |           11           | Kernel size for morphological closing after edge detection to close the remaining contours (only for `method` *canny*) |
|   `min_filter`    |      Flag       |         false          | Use minimum filtering only (only for `method` *canny*)       |
| `min_med_filter`  |      Flag       |         false          | If the image is very noisy, use minimum filtering in addition to the median filter (only for `method` *canny*) |
|   `remove_grid`   |      Flag       |          none          | Method to remove static grid lines from the images (only for `method` *faf*, see **Further remarks**) |
|  `fft_mask_path`  |      Flag       | /data/.../fft_mask.png | Path for fft-mask (only for `method` *faf* and if `remove_grid`=*fft_mask*) |
|  `clear_border`   |      Flag       |         false          | Remove objects that are connected to the edge of the image   |
|  sd1_kernelSize   |   Integer > 0   |           3            | Kernel radius for morphology *Blackhat* operation which is the difference between of closing of the input image and input image (only for ` method` *faf*). For further information see:      https://docs.opencv.org/master/d9/d61/tutorial_py_morphological_ops.html |
|  sd2_kernelSize   |   Integer > 0   |           5            | Kernel radius for morphology *Top Hat* operation which is the difference between of input image and opening of the input image (only for ` method` *faf*). For further information see:   https://docs.opencv.org/master/d9/d61/tutorial_py_morphological_ops.html |
| threshold_binary  | Float = [0,255] |           2.           | Threshold after creation of the standard deviation map, depending on the microscopic technique used, signal and background / noise may be distributed differently. To determine this parameter individually select **`debug`** = `True` (only for ` method` *faf*, see **Further remarks**) |
| erode_kernelSize  |  Integer >= 0   |           3            | Kernel size for the morphology function to erode the final binary segmented objects (only for ` method` *faf*, see **Further remarks**) |
|      `debug`      |      Flag       |         false          | Produces an additional output                                |
|    `n_threads`    |   Integer > 0   |           1            | Number of threads used                                       |

------

## Parameter instructions:

![](./../doc/images/manual_segmentation.png)

------

## Further remarks:

- **`method`** for the segmentation, see the valid types in the following and the config.json file:
  
- `gmm`: Normal PMNs (for brightfield and gray scaled images)
  
- `faf`: Normal PMNs (for brightfield and gray scaled images, method includes flat cells)
  
- `gmmFCs`: Normal PMNs, requires the additional path `input_FCs`  to fungal cells on green fluorescence channel  (for brightfield and gray scaled images, especially for images with flat cells) 
  
  ![](./../doc/images/segmentation_brightfield.png)
  
- `canny`: Fungal cells on green fluorescence channel
  
  ![](./../doc/images/segmentation_green_canny.png)
  
- `singlecells`: Fungal cells on green fluorescence channel (including cluster splitting for convex shape) 
  
  ![](./../doc/images/segmentation_green_singlecells.png)
  
- `otsu`: Cells on red fluorescence channel with dead PI-stained cells
  
    ![](./../doc/images/segmentation_red_otsu.png)
  
- **`remove_grid`** function to remove static grid lines from the images, if necessary (only for method *faf*)
  
- **`remove_grid`** function to remove static grid lines from the images, if necessary (only for method *faf*)

  - `hough` : Usage of hough transformation, which automatically detect grid lines

  - `fft_mask`: remove grid lines in Fourier space by point-wise multiplication of transformed image with a manually created binary mask (preferably created with image processing program Fiji), specify path to binary mask with parameter `fft_mask_path`

    ![](./../doc/images/segmentation_fftMask.png)

- **`threshold_binary`** threshold value after creation of the standard deviation map, depending on the microscopic technique used, signal and background / noise may be distributed differently. To determine this parameter individually select **`debug`** = `True` (only for method *faf*).

    ![](./../doc/images/segmentation_threshold_binary.png)

- **`erode_kernelSize`** kernel size for the morphology function to erode the final binary segmented objects. With large kernel size, small objects can disappear. This functionality can be used to separate connected objects (only for method *faf*).![](./../doc/images/segmentation_erode_kernelSize.png)

