# 2D Wadell roundness and sphericity in Python

Python implementation of Zheng and Hryciw [1] algorithm for Wadell roundness and sphericity calculation. Original implementation in Matlab can be found [here](https://se.mathworks.com/matlabcentral/fileexchange/60651-particle-roundness-and-sphericity-computation)

This version is refactored/updated to enable simpler use and provide more accurate results.

For theoretical details see original paper:
> [1]. Zheng, J., and R. D. Hryciw. “Traditional Soil Particle Sphericity, Roundness and Surface Roughness by Computational Geometry.” Geotechnique, vol. 65, no. 6, 2015, pp. 494–506, https://doi.org/10.1680/geot.14.P.192. 

**For in-depth introduction to the library, visit the [CodeOcean capsule]() (TODO)**

## Installation
```
pip install git+https://github.com/PaPieta/wadell_rs.git
```

## Prerequisites

Library requires ```numpy```, ```scipy```, ```scikit-image```, ```scikit-learn``` and ```edt```. They can be installed via provided requirements file:
```sh
  pip install -r requirements.txt
```


## How to use

The first step is always to load, and pre-process the image.
``` python
from skimage.measure import label
import edt
from wadell_rs import common

img = ... # Load/provide binary image
label_img = label(img) # Assign labels to binary objects
edt_img = edt.edt(img) # Calculate distance transform

obj_dict_list = common.characterize_objects(label_img, edt_img) # Collect characteristics about objects
```

### Sphericity

Choose one of the objects, and sphericity methods and calculate.

``` python
from wadell_rs import sphericity

obj_dict = obj_dict_list[0]
area_sphericity = sphericity.calculate_sphericity(obj_dict, method="area")
```
Available methods:  ```area```, ```diameter```, ```circle_ratio```, ```perimeter```, ```width_to_length```. See documentation or [1] for details.

### Roundness

Choose one of the objects, define parameters and run.

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

For in-depth introduction, visit the [CodeOcean capsule]() (TODO).

## Differences compared to original Matlab implementation

TODO

## External sources

* Sphere fitting to points has been sourced from [sciki-guess library](https://gitlab.com/madphysicist/scikit-guess/-/tree/master)
* Boundary point extraction has been sourced from the [Deepwings project](https://github.com/machine-shop/deepwings/tree/master)
* Smallest enclosing circle calculation has been sourced from [Project Nayuki](https://www.nayuki.io/page/smallest-enclosing-circle)

## License

MIT License (see LICENSE file), excluding external code (see wadell_rs/external for license details in each file).