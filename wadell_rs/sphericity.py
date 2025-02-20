import numpy as np


def calculate_sphericity(object_dict, method="area"):
    """Calculate 2D sphericity of an object.
    Arguments:
        object_dict : dict
            Dictionary with object properties. Output of the function wadell_rs.common.characterize_objects
        method : str
            Method to calculate sphericity. Options are: "area", "diameter", "circle_ratio", "perimeter", "width_to_length"
    Returns:
        float
            2D sphericity value
    """
    if method == "area":
        return area_sphericity(object_dict["area"], object_dict["R_circum"])

    elif method == "diameter":
        return diameter_sphericity(object_dict["area"], object_dict["R_circum"])

    elif method == "circle_ratio":
        return circle_ratio_sphericity(object_dict["R_max"], object_dict["R_circum"])

    elif method == "perimeter":
        return perimeter_sphericity(object_dict["area"], object_dict["perimeter"])

    elif method == "width_to_length":
        return width_to_length_sphericity(object_dict["d1d2"])
    else:
        raise ValueError("Invalid method")


def area_sphericity(object_area, R_circum):
    """2D sphericity defined as a ratio of object area
    to the area of a minimum circumscribing circle.
    Arguments:
        object_area : float
            Area of the object
        R_circum : float
            Radius of the minimum circumscribing circle
    Returns:
        float
            2D sphericity value
    """
    return object_area / (np.pi * R_circum**2)


def diameter_sphericity(object_area, R_circum):
    """2D sphericity defined as a ratio of the diameter of a circle with the same area as the object
    to the diameter of the minimum circumscribing circle.
    Arguments:
        object_area : float
            Area of the object
        R_circum : float
            Radius of the minimum circumscribing circle
    Returns:
        float
            2D sphericity value
    """
    return np.sqrt(object_area / np.pi) / R_circum


def circle_ratio_sphericity(R_max, R_circum):
    """2D sphericity defined as a ratio of the maximum radius of the object
    to the radius of the minimum circumscribing circle.
    Arguments:
        R_max : float
            Maximum radius of the object
        R_circum : float
            Radius of the minimum circumscribing circle
    Returns:
        float
            2D sphericity value
    """
    return R_max / R_circum


def perimeter_sphericity(object_area, object_perimeter):
    """2D sphericity defined as a ratio of the perimeter of a circle with the same area as the object
    to the perimeter of the object. This measure is also known as ISO(2008) circularity.
    Arguments:
        object_area : float
            Area of the object
        object_perimeter : float
            Perimeter of the object
    Returns:
        float
            2D sphericity value
    """
    return 2 * np.sqrt(np.pi * object_area) / object_perimeter


def width_to_length_sphericity(d1d2):
    """2D sphericity defined as a ratio of the with to lehgth of the object.
    Arguments:
        d1d2 : tuple
            Tuple with the object width and length
    Returns:
        float
            2D sphericity value
    """
    return d1d2[1] / d1d2[0]
