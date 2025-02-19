import numpy as np

from skimage.measure import regionprops

from .external import boundary
from .external import smallestenclosingcircle as sec
from . import util


def characterize_objects(label_img, dist_img):
    """Collects properties and boundary points of the objects in the image.
    Called "discrete_boundary" in the original code.
    Arguments:
        label_img: np.ndarray
            Labeled image of the particles
        dist_img: np.ndarray
            Distance transform of the image
    Returns:
        list
            List of dictionaries with properties and boundary points of the objects
            polar - polar coordinates of the boundary points
            rawXY - raw boundary points
            centroid - centroid of the object (see regionprops)
            orientation - orientation of the object (see regionprops)
            area - area of the object (see regionprops)
            perimeter - perimeter of the object (see regionprops)
            d1d2 - major and minor axis lengths of the object (see regionprops major_axis_length and minor_axis_length)
            R_max - radius of the maximum inscribed circle
            R_max_pos - position of the maximum inscribed circle
    """
    props = regionprops(label_img)

    obj_list = []

    for prop in props:
        dist_img_prop = np.copy(dist_img)
        dist_img_prop[label_img != prop.label] = 0

        boundary_vals = boundary.boundary_tracing(prop)

        # Calculate the circumcircle radius
        boundary_vals_tuple = [tuple(row) for row in boundary_vals.tolist()]
        R_circum = sec.make_circle(boundary_vals_tuple)[2]

        polar = util.cart2pol(boundary_vals, prop.centroid)

        obj_dict = {
            "polar": polar,
            "rawXY": np.flip(boundary_vals, axis=1),
            "centroid": prop.centroid,
            "orientation": prop.orientation,
            "area": prop.area,
            "perimeter": prop.perimeter,
            "d1d2": np.array([prop.major_axis_length, prop.minor_axis_length]),
            "R_max": np.max(dist_img_prop),
            "R_max_pos": np.unravel_index(np.argmax(dist_img_prop), dist_img.shape),
            "R_circum": R_circum,
        }

        obj_list.append(obj_dict)

    return obj_list


def curve_smoothing(coords, method="energy", **kwargs):
    """Smoothing of the boundary curve using either loess or energy minimization [1,2]. Named nonparametric_fitting in the original code.
    Original code uses loess smoothing. Energy minimization returns cleaner results but can be slower.
    Arguments:
        coords: (n,2) np.array
            Boundary coordinates
        method: str
            'loess' or 'energy',
                if 'loess' define 'span' and provide polar coordinates in 'coords' and obj_dict['centroid'] as 'centroid' in kwargs,
                if 'energy' define 'alpha_ratio' and 'beta_ratio', and provide rawXY coordinates in 'coords' and obj_dict['perimeter'] as 'perimeter' in kwargs
        kwargs: dict
            Additional arguments for the chosen method
    Returns:
        smoothed: (n,2) np.array
            Smoothed boundary curve
    Examples:
        curve_smoothing(obj_dict['polar'][k], method='loess', span=0.07, centroid=obj_dict['centroid'][k])
        curve_smoothing(obj_dict['rawXY'][k], method='energy', alpha_ratio=0.05, beta_ratio=0.001, perimeter=obj_dict['perimeter'][k])
    References:
    [1]. KASS, M., et al. “SNAKES - ACTIVE CONTOUR MODELS.”
        International Journal of Computer Vision, vol. 1, no. 4, 1987, pp. 321–31,
        https://doi.org/10.1007/BF00133570.
    [2]. Xu, Chenyang & Pham, Dzung & Prince, Jerry. (2000). Image Segmentation Using Deformable Models.
        Handbook of Medical Imaging: Volume 2. Medical Image Processing and Analysis.
    """

    if method == "loess":
        span = kwargs.get("span", 0.07)
        centroid = kwargs.get("centroid", None)

        if centroid is None:
            raise ValueError("Centroid must be provided for loess smoothing")

        coords[:, 0] = util.loess_smoothing(coords[:, 1], coords[:, 0], span)

        smoothed = util.pol2cart(coords, centroid)
        smoothed[-1] = smoothed[0]
        smoothed = np.flip(smoothed, axis=1)

    elif method == "energy":
        alpha_ratio = kwargs.get("alpha_ratio", 0.05)
        beta_ratio = kwargs.get("beta_ratio", 0.001)
        perimeter = kwargs.get("perimeter", None)

        if perimeter is None:
            raise ValueError("Perimeter must be provided for energy minimization")

        N = len(coords) - 1

        alpha = alpha_ratio * perimeter
        beta = beta_ratio * perimeter

        # Create regularization matrix and smooth with backwards Euler
        smoothed = np.matmul(util.regularization_matrix(N, alpha, beta), coords[:-1])

    return smoothed
