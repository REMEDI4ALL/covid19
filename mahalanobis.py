"""
Functions to compute Singular Value Decomposition (SVD), Whitening Tranformation, and Mahalanobis Distance.
Whitening is a transformation that takes the distribution of Healthy cells and maps them to such space that
1. has orthogonal directions – each feature is uncorrelated with others,
2. these directions span the most of the covariance – they capture the scattering of Healthy cells best,
3. new space has lower dimension - we can drop those orthogonal directions that span too little information (same as in PCA),
4. orthogonal directions are standartized – the scattering of Healthy cells becomes uniform along each direction; hence each feauture impacts our distance equally.
Mahalanobis distance = Ordinary Euclidean distance computed after whitening transformation
"""

import numpy as np
from scipy import linalg, stats
# from sklearn.covariance import MinCovDet


def compute_whitening_transform(cells, rank=None):
    cov_matrix = np.cov(cells, rowvar=False)
    sigma, u = linalg.eigh(cov_matrix, check_finite=True)
    descending_idx = sigma.argsort()[::-1]
    sigma = sigma[descending_idx]
    u = u[:, descending_idx]
    
    if rank is None:
        eps = np.finfo(u.dtype.char.lower()).eps
        rtol = np.max(np.abs(sigma)) * max(cov_matrix.shape) * eps
        cutoff = (abs(sigma) > rtol)
        rank = np.sum(cutoff)
    else:
        cutoff = (sigma > 0.)
        rank = min(np.sum(cutoff), rank)
    
    psigma_diag = 1.0 / sigma[:rank]
    whitening_transform = np.diag(np.sqrt(psigma_diag)) @ u.T[:rank, :]
    return whitening_transform, rank


# def compute_whitening_transform_robust(cells, rank=None):
#     u, sigma, vt = linalg.svd(cells - np.mean(cells, axis=0), check_finite=True)
#     descending_idx = sigma.argsort()[::-1]
#     sigma = sigma[descending_idx]
#     u = u[:, descending_idx]
#     vt = vt[descending_idx, :]

#     if rank is None:
#         eps = np.finfo(u.dtype.char.lower()).eps
#         rtol = np.max(np.abs(sigma)) * max(cells.shape) * eps
#         cutoff = (abs(sigma) > rtol)
#         rank = np.sum(cutoff)
#     else:
#         cutoff = (sigma > 0.)
#         rank = min(np.sum(cutoff), rank)
#     cells = u[:, :rank] @ np.diag(sigma[:rank]) @ vt[:rank, :]

#     cov_matrix = MinCovDet().fit(cells).covariance_
#     whitening_transform = linalg.inv(cov_matrix)
#     return whitening_transform, rank


def get_distance(healthy_cells, cells, n_components=None, robust=False):
    if robust:
        # whitening_transform, n_components = compute_whitening_transform_robust(healthy_cells, rank=n_components)
        raise NotImplementedError("To test the variant with robust covariance estimation, you can try the commented-out method above.")
    else:
        whitening_transform, n_components = compute_whitening_transform(healthy_cells, rank=n_components)

    whitening_transform, n_components = compute_whitening_transform(healthy_cells, rank=n_components)
    healthy_mean = np.mean(healthy_cells, axis=0)
    delta = (cells - healthy_mean).T
    delta = whitening_transform @ delta
    distance = np.linalg.norm(delta, ord=2, axis=0)
    return distance


def get_distance_deviation(healthy_cells, cells, n_components=None, robust=False):
    if robust:
        # whitening_transform, n_components = compute_whitening_transform_robust(healthy_cells, rank=n_components)
        raise NotImplementedError("To test the variant with robust covariance estimation, you can try the commented-out method above.")
    else:
        whitening_transform, n_components = compute_whitening_transform(healthy_cells, rank=n_components)

    whitening_transform, n_components = compute_whitening_transform(healthy_cells, rank=n_components)
    healthy_mean = np.mean(healthy_cells, axis=0)
    delta = (cells - healthy_mean).T
    delta = whitening_transform @ delta
    distance = np.linalg.norm(delta, ord=2, axis=0)
    chi_squared_mode = n_components - 2
    deviation = distance - np.sqrt(chi_squared_mode)
    return deviation


def get_proba(healthy_cells, cells, n_components=None, robust=False):
    if robust:
        # whitening_transform, n_components = compute_whitening_transform_robust(healthy_cells, rank=n_components)
        raise NotImplementedError("To test the variant with robust covariance estimation, you can try the commented-out method above.")
    else:
        whitening_transform, n_components = compute_whitening_transform(healthy_cells, rank=n_components)

    whitening_transform, n_components = compute_whitening_transform(healthy_cells, rank=n_components)
    healthy_mean = np.mean(healthy_cells, axis=0)
    delta = (cells - healthy_mean).T
    delta = whitening_transform @ delta
    distance = np.linalg.norm(delta, ord=2, axis=0)
    proba = 1 - stats.chi2.cdf(distance**2, df=n_components)
    return proba


def get_singular_values(X):
    '''Compute singular values and the corresponding rank for observations in X.'''
    cov_matrix = np.cov(X, rowvar=False)
    sigma, u = linalg.eigh(cov_matrix, check_finite=True)
    rank = min(*X.shape) - 1
    sigma = sigma[::-1]
    return sigma, rank



# There are some methods in literature that allow to automatically select the optimal number of components from PCA
# We employ one such method below
def svht(X, sigma=None):
    """
    Return the optimal singular value hard threshold (SVHT) value from [1].
    `X` is any m-by-n matrix. `sigma` is the standard deviation of the 
    noise, if known. Implementation from [2].
    [1] Gavish, Matan, and David L. Donoho. “The optimal hard threshold for singular values is.” IEEE Transactions on Information Theory 60.8 (2014): 5040-5053.
    [2] Robert Taylor. "Optimal Singular Value Hard Threshold." https://humaticlabs.com/blog/optimal-svht/.
    """

    def omega_approx(beta):
        '''Return an approximate omega value for given beta. Equation (5) from Gavish 2014.'''
        return 0.56 * beta**3 - 0.95 * beta**2 + 1.82 * beta + 1.43

    def lambda_star(beta):
        '''Return lambda star for given beta. Equation (11) from Gavish 2014.'''
        return np.sqrt(2 * (beta + 1) + (8 * beta) / 
                    (beta + 1 + np.sqrt(beta**2 + 14 * beta + 1)))
    
    try:
        m, n = sorted(X.shape) # ensures m <= n
    except:
        raise ValueError('invalid input matrix')
    beta = m / n # ratio between 0 and 1
    if sigma is None: # sigma unknown
        sv = linalg.svdvals(X)
        sv = np.squeeze(sv)
        if sv.ndim != 1:
            raise ValueError('vector of singular values must be 1-dimensional')
        return np.median(sv) * omega_approx(beta)
    else: # sigma known
        return lambda_star(beta) * np.sqrt(n) * sigma
