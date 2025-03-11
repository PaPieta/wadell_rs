import numpy as np

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures


def cart2pol(cart, center):
    """Performs a conversion from carthesian to polar coordinates
    Arguments:
        cart:(n,2) np.array
            cartesian coordinates, row/col or y/x format
    Returns: (n,2) np.array
            polar coordinates, (R, phi)
    """

    pol = np.zeros(cart.shape)
    cart = cart - center
    pol[:, 0] = np.sqrt(cart[:, 0] ** 2 + cart[:, 1] ** 2)
    pol[:, 1] = np.arctan2(cart[:, 0], cart[:, 1])
    return pol


def pol2cart(pol, center):
    """Performs a conversion from polar to carthesian coordinates
    Arguments:
        pol:(n,2) np.array
            polar coordinates, (R, phi)
    Returns: (n,2) np.array
            cartesian coordinates, row/col or y/x format
    """

    cart = np.zeros(pol.shape)
    cart[:, 0] = pol[:, 0] * np.sin(pol[:, 1])
    cart[:, 1] = pol[:, 0] * np.cos(pol[:, 1])
    return cart + center


def loess_smoothing(x, y, frac=0.1):
    """Smooth a curve using loess nonparametric fitting. Attempts to mimic Matlab's smooth function.
    Arguments:
        x:(n,) np.array
            x data points
        y:(n,) np.array
            y data points
        frac:float
            fraction of data points to use as neighbourhood for fitting, 0-1 range
    Returns: (n,) np.array
            smoothed y values
    """

    n = len(x)
    smoothed_y = np.zeros(n)
    half_window = int(np.ceil(frac * n / 2))  # Window size based on fraction

    for i in range(n):
        # Find the window around the current point
        left = max(0, i - half_window)
        right = min(n, i + half_window + 1)

        # Select the local data for fitting
        x_window = x[left:right]
        y_window = y[left:right]

        # Fit a 2nd-degree polynomial to the local data
        poly = PolynomialFeatures(degree=2)
        X_poly = poly.fit_transform(x_window.reshape(-1, 1))
        model = LinearRegression().fit(X_poly, y_window)

        # Predict the smoothed value at the current x
        smoothed_y[i] = model.predict(poly.transform([[x[i]]]))[0]

    return smoothed_y


def regularization_matrix(N, alpha, beta):
    """An NxN matrix for imposing elasticity and rigidity to snakes.
    Credits to Vedrana and Anders Dahl course notes https://www2.imm.dtu.dk/courses/02506/
    Arguments:
        N:int
            size of the matrix
        alpha:float
            weight for second derivative (elasticity)
        beta:float
            weight for fourth derivative (rigidity)
    Returns: (N,N) np.array
            regularization matrix
    """

    d = alpha * np.array([-2, 1, 0, 0]) + beta * np.array([-6, 4, -1, 0])
    D = np.fromfunction(
        lambda i, j: np.minimum((i - j) % N, (j - i) % N), (N, N), dtype=int
    )
    A = d[np.minimum(D, len(d) - 1)]
    return np.linalg.inv(np.eye(N) - A)
