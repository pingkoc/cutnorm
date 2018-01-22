import time
import math
import numpy as np
import OptManiMulitBallGBB


def cutnorm_round(U: np.ndarray,
                  V: np.ndarray,
                  C: np.ndarray,
                  max_round_iter: int,
                  logn_lowrank=False,
                  debug=False) -> (np.float_, np.ndarray, np.ndarray, dict):
    '''
    Gaussian Rounding for Cutnorm

    The algorithm picks a random standard multivariate gaussian vector
    w in R^p and computes the rounded solution based on sgn(w \dot ui).

    Adopted from David Koslicki's cutnorm rounding code
    https://github.com/dkoslicki/CutNorm
    and Peter Diao's modifications

    Args:
        U: ndarray, (p, n) shaped matrices of relaxed solutions
        V: ndarray, (p, n) shaped matrices of relaxed solutions
        C: ndarray, original (n, n) shaped matrix to compute cutnorm
        max_round_iter: maximum number of rounding operations
        logn_lowrank: boolean to toggle log2(n) low rank approximation
        debug: boolean, generate debug information
    Returns:
        (approx_opt, uis_opt, vjs_opt, round_debug_info)
        approx_opt: approximated objective function value
        uis_opt: rounded u vector
        vis_opt: rounded v vector
        round_debug_info: debug information for rounding
    '''
    (p, n) = U.shape
    approx_opt = 0
    uis_opt = np.zeros(n)
    vjs_opt = np.zeros(n)
    G = np.random.randn(max_round_iter, p)

    # Debug information
    round_debug_info = {}
    if debug:
        approx_list = np.zeros(max_round_iter)
        uis_list = np.zeros((max_round_iter, n))
        vjs_list = np.zeros((max_round_iter, n))

    # Decomposition for low rank
    if logn_lowrank:
        C_U, C_s, C_V = np.linalg.svd(C)
        low_rank = int(np.log2(n))

    for i in range(max_round_iter):
        g = G[i]
        uis = np.sign(g @ U)
        vjs = np.sign(g @ V)

        # Rounding
        if logn_lowrank:
            C_U_low_filtered = C_U[:, :low_rank] * uis[:, np.newaxis]
            C_V_low_filtered = C_V[:low_rank] * vjs
            C_S = np.diag(C_s[:low_rank])
            C_low_filtered = np.dot(C_U_low_filtered,
                                    np.dot(C_S, C_V_low_filtered))
            approx = np.abs(np.sum(C_low_filtered))
        else:
            approx = np.abs(np.sum(C * np.outer(uis, vjs)))

        if approx > approx_opt:
            approx_opt = approx
            uis_opt = uis
            vjs_opt = vjs

        if debug:
            approx_list[i] = approx / 4.
            uis_list[i] = uis
            vjs_list[i] = vjs

    # Cutnorm is 1/4 of infinity norm
    approx_opt = approx_opt / 4.

    if debug:
        round_debug_info = {
            "round_approx_list": approx_list,
            "round_uis_list": uis_list,
            "round_vjs_list": vjs_list,
            "round_uis_opt": uis_opt,
            "round_vjs_opt": vjs_opt
        }

    return approx_opt, uis_opt, vjs_opt, round_debug_info


def cutnorm_sets(uis: np.ndarray, vjs: np.ndarray) -> (np.ndarray, np.ndarray):
    """
    Generates the cutnorm sets from the rounded SDP solutions

    Args:
        uis: ndarray, (n+1, ) shaped array of rounded +- 1 solution
        vis: ndarray, (n+1, ) shaped array of rounded +- 1 solution
    Returns:
        Reconstructed S and T sets that are {1, 0}^n
        (S, T)
        S: Cutnorm set axis = 0
        T: Cutnorm set axis = 1
    """
    S = -1 * uis[-1] * uis[:-1]
    T = -1 * vjs[-1] * vjs[:-1]

    S = (S + 1) / 2
    T = (T + 1) / 2
    return S, T


def cutnorm(A: np.ndarray,
            B: np.ndarray,
            w1=None,
            w2=None,
            max_round_iter=1000,
            logn_lowrank=False,
            debug=False):
    """
    Computes the cutnorm of the differences between the two matrices

    Args:
        A: ndarray, (n, n) matrix
        B: ndarray, (m, m) matrix
        w1: ndarray, (n, 1) array of weights for A
        w2: ndarray, (m, 1) array of weights for B
        max_round_iter: int, number of iterations for gaussian rounding
        logn_lowrank: boolean to toggle log2(n) low rank approximation
        debug: boolean, generate debug information
    Returns:
        (perf, perf2, S, T, w)
        perf: results from OptManiMulitBallGBB
            [n, p, objf, tsolve, itr, nfe, feasi, nrmG]
            n: dimension of matrix
            p: rank
            objf: objective function value
            tsolve: computation time
            itr, nfe, feasi, nrmG: information from OptManiMulitBallGBB
        perf2: results from gaussian rounding
            [objf_lower, tsolve]
            objf_lower: objective function value from gaussian rounding
            tsolve: computation time
        S: Cutnorm set axis = 0
        T: Cutnorm set axis = 1
        w: weight vector
    Raises:
        ValueError: if A and B are of wrong dimension, or if weight vectors
            does not match the corresponding A and B matrices
    """
    # Input checking
    if A.ndim != 2 or B.ndim != 2:
        raise ValueError("A and B must be 2D matrices")
    n, n2 = np.shape(A)
    m, m2 = np.shape(B)
    if n != n2 or m != m2:
        raise ValueError("A and B must be square matrices")
    if (w1 is None and w2 is not None) or (w1 is not None and w2 is None):
        raise ValueError("Weight vectors required for both matrices")
    if (w1 is not None and w2 is not None and (n != len(w1) or m != len(w2))):
        raise ValueError("Weight vectors need to have the same lenght "
                         "as the first dimension of the corresponding "
                         "matrices")

    if w1 is not None:
        w, C = _compute_C_weighted(A, B, w1, w2)
    else:
        if n == m:
            w, C = _compute_C_eqdim_unweighted(A, B)
        else:
            w, C = _compute_C_uneqdim_unweighted(A, B)

    objf_lower, S, T, debug_info = _compute_cutnorm(C, max_round_iter,
                                                    logn_lowrank, debug)

    return objf_lower, S, T, w, debug_info


def _compute_C_weighted(A: np.ndarray, B: np.ndarray, w1: np.ndarray,
                        w2: np.ndarray) -> (np.ndarray, np.ndarray):
    """
    Generates the difference matrix of the two weighted matricies

    Args:
        A: ndarray, (n, n) matrix
        B: ndarray, (m, m) matrix
        w1: ndarray, (n, 1) array of weights for A
        w2: ndarray, (m, 1) array of weights for B
    Returns:
        C: ndarray, the difference matrix
    """
    v1 = np.hstack((0, np.cumsum(w1)[:-1], 1.))
    v2 = np.hstack((0, np.cumsum(w2)[:-1], 1.))
    v = np.unique(np.hstack((v1, v2)))
    w = np.diff(v)
    n1 = len(w)

    a = np.zeros(n1, dtype=np.int32)
    b = np.zeros(n1, dtype=np.int32)
    for i in range(n1 - 1):
        val = (v[i] + v[i + 1]) / 2
        a[i] = np.argwhere(v1 > val)[0] - 1
        b[i] = np.argwhere(v2 > val)[0] - 1
    # Last element is always itself in new weights
    a[-1] = len(w1) - 1
    b[-1] = len(w2) - 1
    A_sq = A[a]
    A_sq = A_sq[:, a]
    B_sq = B[b]
    B_sq = B_sq[:, b]
    C = (A_sq - B_sq)

    # Normalize C according to weights
    C = C * (np.outer(w, w))
    return w, C


def _compute_C_eqdim_unweighted(A: np.ndarray,
                                B: np.ndarray) -> (np.ndarray, np.ndarray):
    """
    Generates the difference matrix of the two equal dimension unweighted matrices

    Args:
        A: ndarray, (n, n) matrix
        B: ndarray, (m, m) matrix
    Returns:
        C: ndarray, the difference matrix
    """
    n, n2 = np.shape(A)
    n1 = n
    w = np.ones(n) / n
    C = (A - B) / (n1 * n1)  # Normalized C
    return w, C


def _compute_C_uneqdim_unweighted(A: np.ndarray,
                                  B: np.ndarray) -> (np.ndarray, np.ndarray):
    """
    Generates the difference matrix of the two equal dimension unweighted matrices

    Args:
        A: ndarray, (n, n) matrix
        B: ndarray, (m, m) matrix
    Returns:
        C: ndarray, the difference matrix
    """
    n, n2 = np.shape(A)
    m, m2 = np.shape(B)
    d = math.gcd(n, m)
    k = n / d
    l = m / d
    c = k + l - 1
    v1 = np.arange(k) / n
    v2 = np.arange(1, l + 1) / m
    v = np.hstack((v1, v2))
    np.sort(v)
    w = np.diff(v)
    w = np.tile(w, d)

    # Create matrix of differences
    n1 = len(w)
    vals = np.tile(v[:-1], d) + np.floor(np.arange(n1) / c) / d + 1. / (2 * n1)
    a = np.floor(vals * n).astype(int)
    b = np.floor(vals * m).astype(int)
    A_sq = A[a]
    A_sq = A_sq[:, a]
    B_sq = B[b]
    B_sq = B_sq[:, b]
    C = (A_sq - B_sq)

    # Normalize C according to weights
    C = C * (np.outer(w, w))
    return w, C


def _compute_cutnorm(C: np.ndarray,
                     max_round_iter: int,
                     logn_lowrank=False,
                     debug=False) -> (np.float_, np.ndarray, np.ndarray, dict):
    """
    Computes the cutnorm of square matrix C

    Args:
        C: ndarray, (n, n) matrix
        max_round_iter: int, maximum rounding iterations
        logn_lowrank: boolean to toggle log2(n) low rank approximation
        debug: boolean, debug information generation
    Returns:
        (objf_lower, S, T, debug_info)
        objf_lower: objective function value from gaussian rounding
        S: Cutnorm set axis = 0
        T: Cutnorm set axis = 1
        debug_info: dictionary containing debug information
            Debug information from OptManiMulitBallGBB:
                sdp_augm_n: dimension of augmented matrix
                sdp_relax_rank_p: rank
                sdp_objf: objective function value from sdp
                sdp_tsolve: computation time
                sdp_itr, sdp_nfe, sdp_feasi, sdp_nrmG: information from OptManiMulitBallGBB
            Debug information from gaussian rounding:
               round_tsolve: computation time for rounding
               round_approx_list: list of rounded objf values
               round_uis_list: list of uis
               round_vjs_list: list of vjs
               round_uis_opt: optimum uis
               round_vjs_opt: optimum vjs
    """
    n1 = len(C)
    C_col_sum = np.sum(C, axis=0)
    C_row_sum = np.sum(C, axis=1)
    C_tot = np.sum(C_col_sum)
    # Transformation to preserve cutnorm and
    # enforces infinity one norm = 4*cutnorm
    C = np.c_[C, -1.0 * C_row_sum]
    C = np.r_[C, [np.concatenate((-1.0 * C_col_sum, [C_tot]))]]

    # Modify rank estimation for SDP relaxation
    p = int(max(min(round(np.sqrt(2 * n1) / 2), 100), 1))

    # Dim for augmented matrix for SDP
    n2 = 2 * n1 + 2

    # Initial point normalized
    x0 = np.random.randn(p, n2)
    nrmx0 = np.sum(x0 * x0, axis=0)
    x0 = np.divide(x0, np.sqrt(nrmx0))

    tic_sdp = time.time()
    x, g, out = OptManiMulitBallGBB.opt_mani_mulit_ball_gbb(
        x0,
        OptManiMulitBallGBB.cutnorm_quad,
        C,
        record=0,
        mxitr=600,
        gtol=1e-8,
        xtol=1e-8,
        ftol=1e-10,
        tau=1e-3)
    toc_sdp = time.time()
    tsolve_sdp = toc_sdp - tic_sdp

    # SDP upper bound approximation
    U = x[:, :n2 // 2]
    V = x[:, n2 // 2:]
    objf_sdp = np.abs(np.sum(C * (U.T @ V))) / 4.0

    # Gaussian Rounding
    tic_round = time.time()
    (objf_lower, uis, vjs, round_debug_info) = cutnorm_round(
        U, V, C, max_round_iter, logn_lowrank, debug)
    toc_round = time.time()
    tsolve_round = toc_round - tic_round

    # Generate cutnorm sets
    (S, T) = cutnorm_sets(uis, vjs)

    debug_info = {}
    if debug:
        debug_info = {
            "sdp_augm_n": n2,
            "sdp_relax_rank_p": p,
            "sdp_objf": objf_sdp,
            "sdp_tsolve": tsolve_sdp,
            "sdp_itr": out['itr'],
            "sdp_nfe": out['nfe'],
            "sdp_feasi": out['feasi'],
            "sdp_nrmG": out['nrmG'],
            "round_tsolve": tsolve_round
        }
        # Join rounding debug info
        debug_info.update(round_debug_info)

    return objf_lower, S, T, debug_info
