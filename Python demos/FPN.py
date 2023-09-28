import numpy as np
import time
import scipy.linalg
def check_stopping_fpn(opts, x, y, iter, t0):
    # Check the stopping condition of the algorithm
    if 'max_time' not in opts:
        opts['max_time'] = 1e4
    if 'max_iter' not in opts:
        opts['max_iter'] = 1e3
    if 'tol' not in opts:
        opts['tol'] = 1e-6

    stop = 0
    converge = 0
    if time.time() - t0 > opts['max_time']:
        stop = 1  # Maximum CPU time exceeded
    if iter > opts['max_iter']:
        stop = 1  # Maximum iterations exceeded
    if np.linalg.norm(x - y, 'fro') / np.linalg.norm(y, 'fro') < opts['tol']:
        stop = 1  # Condition on successive iterations holds
        converge = 1
    check = {'stop': stop, 'converge': converge}
    return check


def chol_decomposition(X):
    # Check if X is positive definite
    if np.all(np.linalg.eigvals(X) > 0):
        L = np.linalg.cholesky(X)
        flag_pd = 0
    else:
        # Return the partial decomposition like MATLAB
        U, s, V = np.linalg.svd(X)
        rank = np.sum(s > 1e-10)  # adjust the threshold as necessary
        flag_pd = rank + 1  # MATLAB returns the index where decomposition failed
        L = None
    return L, flag_pd


def objective_function(T, S):
    try:
        # Cholesky decomposition
        A = np.linalg.cholesky(T)
        FLAG = 0
    except np.linalg.LinAlgError:
        A = None
        FLAG = 1

    if FLAG == 0:
        lg_det = 2 * np.sum(np.log(np.diag(A)))
        fun = -1 * lg_det + np.dot(T.flatten(), S.flatten())
    else:
        fun = np.inf  # Or any appropriate value indicating the failure

    # Return as a dictionary
    return {"value": fun, "flag": FLAG}

def solver_fpn(S, lmbd, opts=None, hist_num=None):
    t0 = time.time()

    p = S.shape[0]
    iter = 0
    delta = 1e-15
    alpha = 0.5

    if len(lmbd.shape) > 1:
        Lamb = lmbd  # weighted L1 norm
    else:
        Lamb = lmbd * (np.ones((p, p)) - np.eye(p))  # L1 norm

    if hist_num is None:
        hist_num = 0

    opts = opts or {}
    opts.setdefault('max_iter', 1e4)
    opts.setdefault('max_time', 1e4)
    opts.setdefault('tol', 1e-12)
    beta = opts.get('beta', 0.5)
    display = opts.get('display', 1)

    flag_edge = 'edge' in opts
    edgeset = np.where(~np.array(opts.get('edge', []), dtype=bool))[0]

    flag_opt = 'X_opt' in opts
    if flag_opt:
        X_opt = opts['X_opt']
        f_opt = objective_function(X_opt, S - Lamb)['value']
        relobj_iter = []
        relerr_iter = []

    proj = lambda T: np.minimum(0, T - np.diag(np.diag(T))) + np.diag(np.diag(T))

    X = np.diag(1.0 / np.diag(S))
    Fcur = objective_function(X, S - Lamb)
    if Fcur['flag']:
        X = np.diag(1.0 / (np.diag(S) + np.ones(p) * 1e-3))
        Fcur = objective_function(X, S - Lamb)
    objcur = Fcur['value']
    Theta_hist = np.zeros((p, p, 200)) if hist_num > 0 else None
    if Theta_hist is not None and hist_num > 0:
        Theta_hist[:, :, 0] = X
    obj_iter = [objcur]
    time_iter = [time.time() - t0]

    if flag_opt:
        rel_object = abs(f_opt - objcur) / abs(f_opt)
        relobj_iter.append(rel_object)
        rel_err = np.linalg.norm(X_opt - X, 'fro') / np.linalg.norm(X_opt, 'fro')
        relerr_iter.append(rel_err)

    grad = lambda T: -np.linalg.inv(T) + S - Lamb
    gradf = grad(X)

    check = check_stopping_fpn(opts, X, X + 1e16 * np.ones(X.shape), iter, t0)

    if check['stop']:
        X_new = X
        objnew = objcur

    while not check['stop']:
        iter += 1

        rstset = (X - np.diag(np.diag(X)) - 1e3 * np.eye(p) > -delta) & (gradf < 0)

        if flag_edge:
            rstset = np.union1d(rstset, edgeset)

        X_up = X.copy()
        X_up[rstset] = 0

        grady = gradf.copy()
        grady[rstset] = 0

        descent = X_up @ grady @ X_up
        descent[rstset] = 0

        step_size = 1

        Theta_f = lambda gamma: proj(X_up - gamma * descent)
        X_new = Theta_f(step_size)

        A, flag_pd = chol_decomposition(X_new)
        if not flag_pd:
            lg_det = 2 * np.sum(np.log(np.diag(A)))
            objnew = -1 * lg_det + np.dot(X_new.flatten(), (S - Lamb).flatten())
        else:
            objnew = 1e8

        gd = abs(np.dot(gradf.flatten(), descent.flatten()))
        gdI = abs(np.dot(gradf[rstset], X[rstset]))

        # some more computations and updates...

        while (objnew > objcur - alpha * step_size * gd - alpha * gdI) and step_size > np.finfo(float).eps:
            step_size *= beta
            X_new = Theta_f(step_size)

            A, flag_pd = chol_decomposition(X_new)
            if not flag_pd:
                lg_det = 2 * np.sum(np.log(np.diag(A)))
                objnew = -1 * lg_det + np.dot(X_new.flatten(), (S - Lamb).flatten())
            else:
                objnew = 1e8

        check = check_stopping_fpn(opts, X_new, X, iter, t0)

        obj_iter.append(objnew)
        time_iter.append(time.time() - t0)  # Assuming you imported time at the start

        if hist_num > 0 and iter < hist_num:
            # Assuming Theta_hist is pre-initialized and has 3 dimensions
            Theta_hist[:, :, len(obj_iter) - 1] = X_new

        if flag_opt:
            rel_object = abs(f_opt - objnew) / abs(f_opt)
            relobj_iter.append(rel_object)
            rel_err = np.linalg.norm(X_opt - X_new, 'fro') / np.linalg.norm(X_opt, 'fro')
            relerr_iter.append(rel_err)

        if display and iter % 5 == 0:
            print(f"iter: {iter}   objective: {obj_iter[-1]:.11f}   cpu time: {time_iter[-1]:.11f}")

        gradf = grad(X_new)  # Assuming grad is a predefined function
        X = X_new
        objcur = objnew


    run_time = time.time() - t0

    out = {}
    if flag_opt:
        out['relobj_itr'] = relobj_iter
        out['relerr_itr'] = relerr_iter

    out["time"] = run_time
    out["X_est"] = X_new
    out["objective"] = objnew
    out["obj_itr"] = obj_iter
    out["time_itr"] = time_iter
    out["iterate"] = iter
    out["converge"] = check["converge"]

    if hist_num > 0:
        out["Theta_hist"] = Theta_hist[:, :, :len(obj_iter)]
    return out

