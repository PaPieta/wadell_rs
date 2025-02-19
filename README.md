# 2D Wadell roundness and sphericity in Python

Python implementation of Zheng and Hryciw [1] algorithm for Wadell roundness and sphericity calculation. Original implementation in Matlab can be found [here](https://se.mathworks.com/matlabcentral/fileexchange/60651-particle-roundness-and-sphericity-computation).

This version is refactored/updated to enable simpler use and provide more accurate results.

For theoretical details see original paper:
> [1]. Zheng, J., and R. D. Hryciw. “Traditional Soil Particle Sphericity, Roundness and Surface Roughness by Computational Geometry.” Geotechnique, vol. 65, no. 6, 2015, pp. 494–506, https://doi.org/10.1680/geot.14.P.192. 

**For in-depth introduction to the library, visit the [CodeOcean capsule]() (TODO)**

## Installation

``` python
pip install git+https://github.com/PaPieta/wadell_rs.git
```

## Prerequisites

Library requires ```numpy```, ```scipy```, ```scikit-image```, ```scikit-learn``` and ```edt```. They can be installed via provided requirements file:

```sh
  pip install -r requirements.txt
```

## How to use

>For in-depth introduction, visit the [CodeOcean capsule]() (TODO).

The first step is always to load, and pre-process the image.

``` python
from skimage.measure import label
import edt
from wadell_rs import common

img = ... # Load/provide binary image
label_img = label(img) # Assign labels to binary objects
edt_img = edt.edt(img) # Calculate distance transform

obj_dict_list = common.characterize_objects(label_img, edt_img) # Collect characteristics about binary objects
```

### Sphericity

Choose one of the binary objects, and one of the sphericity definitions and calculate.

``` python
from wadell_rs import sphericity

obj_dict = obj_dict_list[0]
area_sphericity = sphericity.calculate_sphericity(obj_dict, method="area")
```

Available definitions:  ```area```, ```diameter```, ```circle_ratio```, ```perimeter```, ```width_to_length```. See documentation or [1] for details.

**Example result:**

<img src="doc_img\sphericity_demo.png" alt="drawing" width="500"/>

### Roundness

Choose one of the binary objects, define parameters and run.

``` python
from wadell_rs import roundness

max_dev_thresh = 0.3 # Maximum deviation from a straight line for discretization
circle_fit_thresh = 0.98 # Defines how close the corner points have to be to the fitted circle outline
smoothing_method = "energy" # Method for smoothing the boundary, 'energy' or 'loess'
alpha_ratio = 0.5 # Ratio of the energy term for the boundary length
beta_ratio = 0.001 # Ratio of the energy term for the boundary curvature

obj_dict = obj_dict_list[0]
roundness_value = roundness.calculate_roundness(
        obj_dict,
        max_dev_thresh,
        circle_fit_thresh,
        smoothing_method=smoothing_method,
        alpha_ratio=alpha_ratio,
        beta_ratio=beta_ratio,
    )
```

**Example result:**

<img src="doc_img\roundness_demo.png" alt="drawing" width="500"/>


## Differences compared to original Matlab implementation

This implementation has a handful of noticeable differences to the original implementation provided in matlab. By far, those cover primarily refactoring of the code, but there is also a few algorithmic differences:

* Originally, the boundary curve was smoothed in ```nonparametric_fitting``` function using Matlab ```smooth(...,'loess')``` function. This approach causes inconsistencies on the ends of the vector creating a boundary, as it doesn't recognize it as a closed loop. Instead, we propose a well-tested image processing method of snake/contour smoothing through energy minimization [2,3]. It can be slower for big boundaries, but seamlessly deals with the closed contours.

* Originally, the ```concave_convex``` function detected middle point position in relation to two edge points and object center. We changed it to simply calculating the angle between three points. This lets us cleanly calculate convexity in one go. No decrease in result quality was observed due to this change

* In circle fitting to corners, there was an issue where sometimes first and last convex points in a boundary were mistakenly assigned two separate circles, and had to be treated separately at the end of the fitting algorithm. We mitigate this problem by re-indexing the convex points so that they start at:
  * A convex point right after a concave section, or
  * The most flat convex section if no concave section is present

* In ```discretize_boundary``` we make sure to correctly handle boundary ends, so that the last two points are not very close together

* Most of the functions are adjusted to run on single objects, and not on all image objects as one. This allows the user to choose objects of interest in case the image is very big or contains noise. The only exception is ```common.characterize_objects```, there we follow the style of methods similar to ```regionprops```, characterizing all objects at once. This part is fast enough to allow that.

>[2]. KASS, M., et al. “SNAKES - ACTIVE CONTOUR MODELS.” International Journal of Computer Vision, vol. 1, no. 4, 1987, pp. 321–31,https://doi.org/10.1007/BF00133570.
>
>[3]. Xu, Chenyang & Pham, Dzung & Prince, Jerry. (2000). Image Segmentation Using Deformable Models.Handbook of Medical Imaging: Volume 2. Medical Image Processing and Analysis.

## External sources

* Sphere fitting to points has been sourced from [scikit-guess library](https://gitlab.com/madphysicist/scikit-guess/-/tree/master)
* Boundary point extraction has been sourced from the [Deepwings project](https://github.com/machine-shop/deepwings/tree/master)
* Smallest enclosing circle calculation has been sourced from [Project Nayuki](https://www.nayuki.io/page/smallest-enclosing-circle)

## License

MIT License (see LICENSE file), excluding external code (see wadell_rs/external for license details in each file).