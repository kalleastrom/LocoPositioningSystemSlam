"""
Functions used to localize the anchors, refer to the publication below
and use the bib-tex entry for citation

@inbook{
    6505724206a7482c8bec0010f5a11197,
    title = "Robust Time-of-Arrival Self Calibration with Missing Data and Outliers",
    author = "Batstone, {Kenneth John} and Magnus Oskarsson and Karl \AA str\o" m",
    year = "2016",
    month = "9",
    doi = "10.1109/EUSIPCO.2016.7760673",
    pages = "2370--2374",
    booktitle = "2016 24th European Signal Processing Conference (EUSIPCO)",
    publisher = "Institute of Electrical and Electronics Engineers Inc.",
    address = "United States"
}
"""
import sys
from math import sqrt
import numpy as np
from scipy import linalg
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsqr
from numpy.linalg import svd, cholesky, inv, eig, solve, norm, lstsq, qr
from numpy.matlib import repmat

def pre_process_compaction_matrix(Bhat, dim):
    u, s, v = svd(Bhat[1:,1:])
    s = np.diag(s)
    xr  = np.transpose(u[:,0:dim])
    yr  = np.dot(s[0:dim,0:dim], np.transpose(v[:,0:dim]))
    xtp = np.concatenate((np.zeros((dim,1)), xr),1)
    yt  = np.concatenate((np.zeros((dim,1)), yr),1)
    xt  = xtp/-2.0;
    return xt, yt

def get_linear_constraints(Bhat, xt):
    # Construct linear constraints
    D = 3
    Nunknowns = int((D * ( D + 1.0 ) ) / 2.0 + D + 1.0)
    cfm_linear = np.zeros((6 - 1, Nunknowns))
    cfm_linear_normalized = np.zeros((6 - 1, Nunknowns))
    
    # Linear equations formed by (-2*xt'*bv + xt'*Cv*xt) - Bhat
    for ii in range(5):
        xr1, xr2, xr3 = xt[:,ii+1]
        cfm_linear[ii,0] = xr1 ** 2.0
        cfm_linear[ii,1] = 2.0 * xr1 * xr2
        cfm_linear[ii,2] = 2.0 * xr1 * xr3
        cfm_linear[ii,3] = xr2 ** 2.0
        cfm_linear[ii,4] = 2.0 * xr2 * xr3
        cfm_linear[ii,5] = xr3 ** 2.0
        cfm_linear[ii,6] = -2.0 * xr1
        cfm_linear[ii,7] = -2.0 * xr2
        cfm_linear[ii,8] = -2.0 * xr3
        cfm_linear[ii,9] = -Bhat[ii+1, 0]

        # normalize constraints
        cfm_linear_normalized[ii,:] = cfm_linear[ii,:]/norm(cfm_linear[ii,:])

    # compute constraints
    AA = cfm_linear_normalized[:,:-1];
    bb = cfm_linear_normalized[:,-1:];
    zz = -linalg.pinv(AA).dot(bb);

    H = np.array([
        [zz[0][0], zz[1][0], zz[2][0]],
        [zz[1][0], zz[3][0], zz[4][0]],
        [zz[2][0], zz[4][0], zz[5][0]]
    ])
    b = np.array([
        [zz[6][0]],
        [zz[7][0]],
        [zz[8][0]]
    ])
    return H, b, cfm_linear

def safe_cholesky_factorization(H):
    realEig = np.real(eig(H)[0])
    if min(realEig) > 0:
        L = cholesky(inv(H))
    else:
        H = H + (-min(realEig) + 0.1) * np.eye(H.shape[0])
        L = cholesky(inv(H))
    return H, L

def toa_3D_bundle(x, y, solution, settings):
    (I, J) = solution.inliers.nonzero()
    D = np.zeros(len(I),)
    #ind = np.ravel_multi_index((I, J), dims=d.shape)
    for n in range(len(I)):
        D[n] = solution.d[I[n],J[n]]
    I = np.reshape(I, (len(I),1))
    J = np.reshape(J, (len(J),1))
    D = np.reshape(D, (len(D),1))
    return bundletoa(D, I, J, x, y, settings)

def calcresandjac(D, I, J, x, y):

    m, n = x.shape[1], y.shape[1]
    V = x[:, I.T[0]] - y[:, J.T[0]]
    Vt = V.T
    dd = np.sqrt([((V ** 2).sum(axis=0))]).T
    idd = 1.0 / dd
    res = dd - D
    
    II = np.array([np.arange(len(I))]).T
    JJ_i = I * 3
    JJ_j = J * 3
    JJ1 = JJ_i + 0
    JJ2 = JJ_i + 1
    JJ3 = JJ_i + 2
    JJ4 = JJ_j + 3 * m
    JJ5 = JJ_j + 3 * m + 1
    JJ6 = JJ_j + 3 * m + 2

    Vt0 = np.array([Vt[:, 0]]).T
    Vt1 = np.array([Vt[:, 1]]).T
    Vt2 = np.array([Vt[:, 2]]).T
    VV1 = idd * Vt0
    VV2 = idd * Vt1
    VV3 = idd * Vt2

    row_ind = np.concatenate((II, II, II, II, II, II)).T[0]
    col_ind = np.concatenate((JJ1, JJ2, JJ3, JJ4, JJ5, JJ6)).T[0]
    data = np.concatenate((VV1, VV2, VV3, -VV1, -VV2, -VV3)).T[0]
    M = len(D)
    N = 3 * m + 3 * n
    jac = csr_matrix((data, (row_ind, col_ind)), shape=(M, N))
    return res, jac

def updatexy(x, y, dz):
    m, n = x.shape[1], y.shape[1]
    dzx = dz[0:(3 * m)]
    dzy = dz[(3 * m):]
    xnew = x + np.reshape(dzx, (m,3)).T
    ynew = y + np.reshape(dzy, (n,3)).T
    return xnew, ynew

def bundletoa(D, I, J, xt, yt, settings):
    
    numIter = settings.bundle.numberOfIterations
    countLim = settings.bundle.counterLimit
    numLim = settings.bundle.numericalLimit

    for iteration in range(numIter):
        sys.stdout.write('\rBundle LS iteration %d' % iteration)
        sys.stdout.flush()
        res, jac = calcresandjac(D, I, J, xt, yt)
        dz_A = -(jac.T.dot(jac) + np.eye(jac.shape[1]))
        dz_b = jac.T.dot(res)
        # Suporior to numpy's in terms of speed and to the sparse solver lsqr
        # in terms of accuracy
        dz = linalg.solve(dz_A, dz_b)
        xtn, ytn = updatexy(xt, yt, dz)
        res2, jac2 = calcresandjac(D, I, J, xtn, ytn)
        
        cc = norm(jac.dot(dz)) / norm(res)
        if norm(res) < norm(res2):
            if cc > numLim:
                counter = 0
                while (counter < countLim) and (norm(res) < norm(res2)):
                    dz = dz / float(2)
                    xtn, ytn = updatexy(xt, yt, dz)
                    res2, _ = calcresandjac(D, I, J, xtn, ytn)
                    counter = counter + 1
        else:
            xt = xtn
            yt = ytn

    return xt, yt, res, jac
    
def toa_3d_bundle(d, x, y, inliers, settings):
    # Extract indices and feasible meaasurements in row vector form
    I, J = np.where(inliers==1)
    D = d[I, J]

    # Convert to make the vecors two-dimensional
    I = np.array([I]).T
    J = np.array([J]).T
    D = np.array([D]).T
        
    # Return solution
    return bundletoa(D, I, J, x, y, settings)

def toa_normalise(x0, y0):
    xdim = x0.shape[0]
    m = x0.shape[1]
    n = y0.shape[1]
    
    # translation
    t = -np.array([x0[:,0]]).T
    x = x0 + repmat(t, 1, m)
    y = y0 + repmat(t, 1, n)

    # rotation
    [q, r]=qr(x[:,1:xdim+1])
    x = np.dot(q.T, x)
    y = np.dot(q.T, y)

    # mirroring
    M = np.diag(np.sign(np.diag(x[:,1:xdim+1])));
    x1 = np.dot(M, x);
    y1 = np.dot(M, y);
    return x1, y1
