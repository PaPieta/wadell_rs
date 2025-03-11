import numpy as np

from .external import nsphere
from . import common


def calculate_roundness(
    obj_dict,
    max_dev_thresh,
    circle_fit_thresh,
    verbose=False,
    smoothing_method="energy",
    **kwargs,
):
    """Calculate the roundness of an object using the Zheng and Hryciw algorithm [1].
    Arguments:
        obj_dict: dict
            Dictionary with properties and boundary points of the object, output of common.characterize_objects
        max_dev_thresh: float
            Maximum deviation from a straight line for discretization
        circle_fit_thresh: float
            Threshold for the circle fitting, defines how close the points have to be to the circle outline, (0-1) range
        verbose: bool
            If True, return additional information in a dict, useful for plotting
        smoothing_method: str
            Method for smoothing the boundary, 'energy' or 'loess'
        **kwargs:
            If smoothing_method is 'energy':
                alpha_ratio: float
                    Ratio of the energy term for the boundary length
                beta_ratio: float
                    Ratio of the energy term for the boundary curvature
            If smoothing_method is 'loess':
                span: float
                    Span for the loess smoothing
    Returns:
        roundess_val: float
            Roundness value
        (if verbose==True) out_dict: dict
            Dictionary with additional information
                R_max: float
                    Radius of the maximum inscribed circle
                R_max_pos: (2,) float np.array
                    Position of the maximum inscribed circle
                boundary_points: (n,2) float np.array
                    Smoothed boundary points
                concave_points: (n,2) float np.array
                    Concave points
                convex_points: (n,2) float np.array
                    Convex points
                corner_circle_radii: (n,) float np.array
                    Radii of the detected corner circles
                corner_circle_centers: (n,2) float np.array
                    Centers of the detected corner circles
                point_idx_groups: list
                    List of used point indices for each corner circle
    References:
    [1] Zheng, J., and R. D. Hryciw.
        “Traditional Soil Particle Sphericity, Roundness and Surface Roughness by Computational Geometry.”
        Geotechnique, vol. 65, no. 6, 2015, pp. 494–506, https://doi.org/10.1680/geot.14.P.192.
    """
    # Smooth the boundary
    if smoothing_method == "energy":
        alpha_ratio = kwargs.get("alpha_ratio", 0.05)
        beta_ratio = kwargs.get("beta_ratio", 0.001)

        boundary_points = common.curve_smoothing(
            obj_dict["rawXY"],
            method="energy",
            alpha_ratio=alpha_ratio,
            beta_ratio=beta_ratio,
            perimeter=obj_dict["perimeter"],
        )

    elif smoothing_method == "loess":
        span = kwargs.get("span", 0.07)

        boundary_points = common.curve_smoothing(
            obj_dict["polar"], method="loess", span=span, centroid=obj_dict["centroid"]
        )

    else:
        raise ValueError('Invalid smoothing method, choose either "energy" or "loess"')

    # Discretize the boundary
    keypoints = discretize_boundary(boundary_points, max_dev_thresh)
    # Find concave and convex points
    concave_points, convex_points = concave_convex(keypoints)
    # Fit circles to the convex points
    radii, centers, point_idx_groups = compute_corner_circles(
        convex_points,
        keypoints,
        obj_dict["R_max"],
        obj_dict["R_max_pos"],
        circle_fit_thresh,
    )

    roundess_val = np.mean(radii) / obj_dict["R_max"]

    if not verbose:
        return roundess_val
    else:
        out_dict = {
            "R_max": obj_dict["R_max"],
            "R_max_pos": obj_dict["R_max_pos"],
            "boundary_points": boundary_points,
            "concave_points": concave_points,
            "convex_points": convex_points,
            "corner_circle_radii": radii,
            "corner_circle_centers": centers,
            "point_idx_groups": point_idx_groups,
        }
        return roundess_val, out_dict


def maxlinedev(X, Y):
    """Find the point with the maximum deviation from a line connecting the first and last points.
    Arguments:
        X:(n,) np.array
            x data points
        Y:(n,) np.array
            y data points
    Returns:
        max_dev: float
            Maximum deviation from the line
        pos_idx: int
            Index of the point with the maximum deviation
    """
    eps = 1e-6

    if len(X) == 1:
        print("Warning: Contour with only one point")
        max_dev = 0
        pos_idx = 0
        return max_dev, pos_idx
    elif len(X) == 0:
        raise ("Error: Contour with no points")

    # End point distance
    dist_end = np.sqrt((X[0] - X[-1]) ** 2 + (Y[0] - Y[-1]) ** 2)

    # If end points are coincidental
    if dist_end < eps:
        # Distance from first point
        dist = np.sqrt((X - X[0]) ** 2 + (Y - Y[0]) ** 2)
    else:
        # Distance from line
        dist = (
            np.abs(
                (Y[0] - Y[-1]) * X + (X[-1] - X[0]) * Y + Y[-1] * X[0] - Y[0] * X[-1]
            )
            / dist_end
        )

    # print(dist)
    max_dev = np.max(dist)
    pos_idx = np.argmax(dist)

    return max_dev, pos_idx


def discretize_boundary(boundary, max_dev_thresh):
    """Discretize the boundary into keypoints based on the maximum deviation from a straight line.
    Arguments:
        boundary: (n,2) float np.array
            Boundary points
        max_dev_thresh: float
            Maximum deviation from a straight line
    Returns:
        keypoints: (m,2) float np.array
            Keypoints filtered from the boundary
    """

    X = boundary[:, 1]
    Y = boundary[:, 0]

    segment_start_idx = 0
    segment_end_idx = len(X)

    keypoints = []
    keypoints.append([X[segment_start_idx], Y[segment_start_idx]])

    # Append the first point to the end of X and Y (ensures cleaner closing of the boundary)
    # (notice that segment_end_idx is not updated for the first loop)
    X = np.append(X, X[0])
    Y = np.append(Y, Y[0])

    while segment_start_idx < segment_end_idx:
        # Find the point with the maximum deviation
        [max_dev, pos_idx] = maxlinedev(
            X[segment_start_idx:segment_end_idx], Y[segment_start_idx:segment_end_idx]
        )

        while max_dev > max_dev_thresh:
            # Update the segment end index
            segment_end_idx = pos_idx + segment_start_idx

            # Find the point with the maximum deviation
            [max_dev, pos_idx] = maxlinedev(
                X[segment_start_idx:segment_end_idx],
                Y[segment_start_idx:segment_end_idx],
            )

        # If we are not at the start, and the last point is the end, skip it as it is already in the list
        if segment_end_idx != len(X) or segment_start_idx == 0:
            # Add the keypoint to the list
            keypoints.append([X[segment_end_idx - 1], Y[segment_end_idx - 1]])
            # print(segment_end_idx)

        # Update the segment indices
        segment_start_idx = segment_end_idx
        segment_end_idx = len(X)

    return np.array(keypoints)


def concave_convex(keypoints):
    """Find the concave and convex points from a list of keypoints.
    Compared to original implementation, this method compares point vector angles, instead of distance to object center.
    Also shifts the keypoints to start from the first concave point to avoid circle fitting issues with the first point.
    Arguments:
        keypoints: (m,2) float np.array
            Keypoints of the boundary
    Returns:
        concave_points: (n,2) float np.array
            Concave points
        convex_points: (n,2) float np.array
            Convex points
    """
    concave_points = []
    convex_points = []

    # Append last point to the beginning
    keypoints = np.vstack([keypoints[-1], keypoints, keypoints[0]])

    v1 = keypoints[1:-1] - keypoints[:-2]
    v2 = keypoints[2:] - keypoints[1:-1]

    angle = np.arctan2(
        v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0],
        v1[:, 0] * v2[:, 0] + v1[:, 1] * v2[:, 1],
    )

    # Find index of the first concave point
    concave_points_loc = np.where(angle > 0)[0]
    if len(concave_points_loc) != 0:
        roll_index = concave_points_loc[0]
    else:
        # Roll to the bigest angle instead
        roll_index = np.argmax(angle)
        # roll_index = 0

    # Take only actual keypoints before rolling
    keypoints = keypoints[1:-1]
    # Shift the keypoints to start from the first concave point
    keypoints = np.roll(keypoints, -roll_index, axis=0)
    angle = np.roll(angle, -roll_index)

    concave_points = keypoints[angle > 0]
    convex_points = keypoints[angle < 0]

    return np.array(concave_points), np.array(convex_points)


def compute_corner_circles(
    convex_points, keypoints, rad_max, rad_max_pos, circle_fit_thresh
):
    """Fits circles to the convex points of the boundary.
    Arguments:
        convex_points: (n,2) float np.array
            Convex points of the boundary
        keypoints: (m,2) float np.array
            Keypoints of the boundary
        rad_max: float
            Radius of the max circle inscribed in the object
        rad_max_pos: (2,) float np.array
            Position of the max circle inscribed in the object
        circle_fit_thresh: float
            Threshold for the circle fitting, defines how close the points have to be to the circle outline, (0-1) range
    Returns:
        radii: (n,) float np.array
            Radii of the detected corner circles
        centers: (n,2) float np.array
            Centers of the detected corner circles
        point_idx_groups: list
            List of used point indices for each corner circle"""

    point_idx_start = 0
    finish_idx = len(convex_points)
    radii = []
    centers = []
    point_idx_groups = []

    # Loop over starting points
    while point_idx_start < finish_idx - 1:
        # Backwards loop over end points
        for point_idx_end in range(finish_idx, point_idx_start + 2, -1):
            # Fit a circle to the points
            radius, center = nsphere.nsphere_fit(
                convex_points[point_idx_start:point_idx_end], scaling=True
            )

            # Distance from the center to the max circle
            dist_circle = np.linalg.norm(center - rad_max_pos)
            # Distance between the boundary and the max circle
            dist_bbox = np.linalg.norm(
                convex_points[point_idx_start:point_idx_end] - rad_max_pos, axis=1
            ).mean()
            if dist_circle < dist_bbox and radius < rad_max:
                # Get distance from circle center to ALL points
                dist = np.linalg.norm(keypoints - center, axis=1)
                # If no point is closer to the circle center than the radius
                # (Within a threshold, fitting is not perfect)
                if np.all(dist / radius > circle_fit_thresh):
                    # Save the circle
                    radii.append(radius)
                    centers.append(center)
                    point_idx_groups.append(np.arange(point_idx_start, point_idx_end))
                    # Update the starting point index
                    point_idx_start = point_idx_end - 1
                    break
        # Update the starting point index
        point_idx_start += 1

    return np.array(radii), np.array(centers), point_idx_groups
