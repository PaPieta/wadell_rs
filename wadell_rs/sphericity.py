import numpy as np


def calculate_sphericity(object_dict, method="area"):
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
    return object_area / (np.pi * R_circum**2)


def diameter_sphericity(object_area, R_circum):
    return np.sqrt(object_area / np.pi) / R_circum


def circle_ratio_sphericity(R_max, R_circum):
    return R_max / R_circum


def perimeter_sphericity(object_area, object_perimeter):
    return 2 * np.sqrt(np.pi * object_area) / object_perimeter


def width_to_length_sphericity(d1d2):
    return d1d2[1] / d1d2[0]
