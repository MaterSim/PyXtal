"""
A minimal collection of optimization algorithms adapted from scipy 
"""

from __future__ import division, print_function, absolute_import
import warnings
import sys
import numpy
from scipy._lib.six import callable, xrange
from numpy import (
    atleast_1d,
    eye,
    mgrid,
    argmin,
    zeros,
    shape,
    squeeze,
    vectorize,
    asarray,
    sqrt,
    Inf,
    asfarray,
    isinf,
)
import numpy as np
from scipy.optimize.linesearch import (
    line_search_wolfe1,
    line_search_wolfe2,
    line_search_wolfe2 as line_search,
    LineSearchWarning,
)

_epsilon = sqrt(numpy.finfo(float).eps)

# standard status messages of optimizers
_status_message = {
    "success": "Optimization terminated successfully.",
    "maxfev": "Maximum number of function evaluations has " "been exceeded.",
    "maxiter": "Maximum number of iterations has been " "exceeded.",
    "pr_loss": "Desired error not necessarily achieved due " "to precision loss.",
}


class _LineSearchError(RuntimeError):
    pass


def _line_search_wolfe12(f, fprime, xk, pk, gfk, old_fval, old_old_fval, **kwargs):
    """
    Same as line_search_wolfe1, but fall back to line_search_wolfe2 if
    suitable step length is not found, and raise an exception if a
    suitable step length is not found.
    Raises
    ------
    _LineSearchError
        If no suitable step size is found
    """

    extra_condition = kwargs.pop("extra_condition", None)

    ret = line_search_wolfe1(f, fprime, xk, pk, gfk, old_fval, old_old_fval, **kwargs)

    if ret[0] is not None and extra_condition is not None:
        xp1 = xk + ret[0] * pk
        if not extra_condition(ret[0], xp1, ret[3], ret[5]):
            # Reject step if extra_condition fails
            ret = (None,)

    if ret[0] is None:
        # line search failed: try different one.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", LineSearchWarning)
            kwargs2 = {}
            for key in ("c1", "c2", "amax"):
                if key in kwargs:
                    kwargs2[key] = kwargs[key]
            ret = line_search_wolfe2(
                f,
                fprime,
                xk,
                pk,
                gfk,
                old_fval,
                old_old_fval,
                extra_condition=extra_condition,
                **kwargs2
            )

    if ret[0] is None:
        raise _LineSearchError()

    return ret


def vecnorm(x, ord=2):
    if ord == Inf:
        return numpy.amax(numpy.abs(x))
    elif ord == -Inf:
        return numpy.amin(numpy.abs(x))
    else:
        return numpy.sum(numpy.abs(x) ** ord, axis=0) ** (1.0 / ord)


def wrap_function(function, args):
    ncalls = [0]
    # print(wrapper_args)
    # sys.exit()
    def function_wrapper(*wrapper_args):
        ncalls[0] += 1
        return function(*(wrapper_args + args))

    return ncalls, function_wrapper


class OptimizeResult(dict):
    """ Represents the optimization result.
    Attributes
    ----------
    x : ndarray
        The solution of the optimization.
    success : bool
        Whether or not the optimizer exited successfully.
    status : int
        Termination status of the optimizer. Its value depends on the
        underlying solver. Refer to `message` for details.
    message : str
        Description of the cause of the termination.
    fun, jac, hess: ndarray
        Values of objective function, its Jacobian and its Hessian (if
        available). The Hessians may be approximations, see the documentation
        of the function in question.
    hess_inv : object
        Inverse of the objective function's Hessian; may be an approximation.
        Not available for all solvers. The type of this attribute may be
        either np.ndarray or scipy.sparse.linalg.LinearOperator.
    nfev, njev, nhev : int
        Number of evaluations of the objective functions and of its
        Jacobian and Hessian.
    nit : int
        Number of iterations performed by the optimizer.
    maxcv : float
        The maximum constraint violation.
    Notes
    -----
    There may be additional attributes not listed above depending of the
    specific solver. Since this class is essentially a subclass of dict
    with attribute accessors, one can see which attributes are available
    using the `keys()` method.
    """

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return "\n".join(
                [k.rjust(m) + ": " + repr(v) for k, v in sorted(self.items())]
            )
        else:
            return self.__class__.__name__ + "()"

    def __dir__(self):
        return list(self.keys())


def _minimize_tpgd(fun, x0, args, jac, beta=1.001, gtol=1e-3, norm=Inf, maxiter=None):
    """
    Gradient descent based on Barzilai-Borwein[1] method
    https://en.wikipedia.org/wiki/Gradient_descent
    """
    f = fun
    fprime = jac
    if maxiter is None:
        maxiter = len(x0) * 200
    k = 0
    func_calls, f = wrap_function(f, args)
    grad_calls, myfprime = wrap_function(fprime, args)
    fx0, gx0 = f(x0), myfprime(x0)
    gnorm = vecnorm(gx0, ord=norm)
    gamma0 = 1e-3 / gnorm
    x1 = x0 + gamma0 * gx0
    fx1, gx1 = f(x1), myfprime(x1)
    gnorm = vecnorm(gx1, ord=norm)

    while (gnorm > gtol) and (k < maxiter):
        (dim, mu, shift) = args
        args = (dim, mu * beta, shift)
        # compute gamma based on the Barzilai-Borwein[1] method
        # https://en.wikipedia.org/wiki/Gradient_descent
        dg = gx1 - gx0
        dx = x1 - x0
        gamma = np.abs(np.dot(dx, dg) / np.dot(dg, dg))

        # update the x0, gx0, fx0
        x0, gx0, fx0 = x1.copy(), gx1.copy(), fx1.copy()
        x1 = x1 - gx1 * gamma
        fx1, gx1 = f(x1), myfprime(x1)
        gnorm = vecnorm(gx1, ord=norm)
        k += 1
        # print('step in tpgd: {:4d} {:12.4f} {:12.4f}'.format(count, fx1, np.max(gx1)))

    result = OptimizeResult(fun=fx1, jac=gx1, nfev=func_calls[0], x=x1, nit=k)
    return result


def _minimize_bfgs(
    fun,
    x0,
    args=(),
    jac=None,
    beta=1.001,
    callback=None,
    gtol=1e-5,
    norm=Inf,
    eps=_epsilon,
    maxiter=None,
    disp=False,
    return_all=False,
    **unknown_options
):
    """
    Minimization of scalar function of one or more variables using the
    BFGS algorithm.
    Options
    -------
    disp : bool
        Set to True to print convergence messages.
    maxiter : int
        Maximum number of iterations to perform.
    gtol : float
        Gradient norm must be less than `gtol` before successful
        termination.
    norm : float
        Order of norm (Inf is max, -Inf is min).
    eps : float or ndarray
        If `jac` is approximated, use this value for the step size.
    """
    f = fun
    fprime = jac
    epsilon = eps
    retall = return_all

    x0 = asarray(x0).flatten()
    if x0.ndim == 0:
        x0.shape = (1,)
    if maxiter is None:
        maxiter = len(x0) * 200
    func_calls, f = wrap_function(f, args)
    grad_calls, myfprime = wrap_function(fprime, args)
    gfk = myfprime(x0)
    k = 0
    N = len(x0)
    I = numpy.eye(N, dtype=int)
    Hk = I

    # Sets the initial step guess to dx ~ 1
    old_fval = f(x0)
    old_old_fval = old_fval + np.linalg.norm(gfk) / 2

    xk = x0
    if retall:
        allvecs = [x0]
    warnflag = 0
    gnorm = vecnorm(gfk, ord=norm)
    while (gnorm > gtol) and (k < maxiter):
        (dim, mu, shift) = args
        args = (dim, mu * beta, shift)

        pk = -numpy.dot(Hk, gfk)
        try:
            alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = _line_search_wolfe12(
                f,
                myfprime,
                xk,
                pk,
                gfk,
                old_fval,
                old_old_fval,
                amin=1e-100,
                amax=1e100,
            )
        except _LineSearchError:
            # Line search failed to find a better solution.
            warnflag = 2
            break

        xkp1 = xk + alpha_k * pk
        if retall:
            allvecs.append(xkp1)
        sk = xkp1 - xk
        xk = xkp1
        if gfkp1 is None:
            gfkp1 = myfprime(xkp1)

        yk = gfkp1 - gfk
        gfk = gfkp1
        if callback is not None:
            callback(xk)
        k += 1
        gnorm = vecnorm(gfk, ord=norm)
        if gnorm <= gtol:
            break

        if not numpy.isfinite(old_fval):
            # We correctly found +-Inf as optimal value, or something went
            # wrong.
            warnflag = 2
            break

        try:  # this was handled in numeric, let it remaines for more safety
            rhok = 1.0 / (numpy.dot(yk, sk))
        except ZeroDivisionError:
            rhok = 1000.0
            if disp:
                print("Divide-by-zero encountered: rhok assumed large")
        if isinf(rhok):  # this is patch for numpy
            rhok = 1000.0
            if disp:
                print("Divide-by-zero encountered: rhok assumed large")
        A1 = I - sk[:, numpy.newaxis] * yk[numpy.newaxis, :] * rhok
        A2 = I - yk[:, numpy.newaxis] * sk[numpy.newaxis, :] * rhok
        Hk = numpy.dot(A1, numpy.dot(Hk, A2)) + (
            rhok * sk[:, numpy.newaxis] * sk[numpy.newaxis, :]
        )

    fval = old_fval
    if np.isnan(fval):
        # This can happen if the first call to f returned NaN;
        # the loop is then never entered.
        warnflag = 2

    if warnflag == 2:
        msg = _status_message["pr_loss"]
    elif k >= maxiter:
        warnflag = 1
        msg = _status_message["maxiter"]
    else:
        msg = _status_message["success"]

    if disp:
        print("%s%s" % ("Warning: " if warnflag != 0 else "", msg))
        print("         Current function value: %f" % fval)
        print("         Iterations: %d" % k)
        print("         Function evaluations: %d" % func_calls[0])
        print("         Gradient evaluations: %d" % grad_calls[0])

    result = OptimizeResult(
        fun=fval,
        jac=gfk,
        hess_inv=Hk,
        nfev=func_calls[0],
        njev=grad_calls[0],
        status=warnflag,
        success=(warnflag == 0),
        message=msg,
        x=xk,
        nit=k,
    )
    if retall:
        result["allvecs"] = allvecs
    return result


def _minimize_cg(
    fun,
    x0,
    args=(),
    jac=None,
    beta=1.001,
    callback=None,
    gtol=1e-5,
    norm=Inf,
    eps=_epsilon,
    maxiter=None,
    disp=False,
    return_all=False,
    **unknown_options
):
    """
    Minimization of scalar function of one or more variables using the
    conjugate gradient algorithm.
    Options
    -------
    disp : bool
        Set to True to print convergence messages.
    maxiter : int
        Maximum number of iterations to perform.
    gtol : float
        Gradient norm must be less than `gtol` before successful
        termination.
    norm : float
        Order of norm (Inf is max, -Inf is min).
    eps : float or ndarray
        If `jac` is approximated, use this value for the step size.
    """
    f = fun
    fprime = jac
    epsilon = eps
    retall = return_all

    x0 = asarray(x0).flatten()
    if maxiter is None:
        maxiter = len(x0) * 200
    func_calls, f = wrap_function(f, args)
    grad_calls, myfprime = wrap_function(fprime, args)

    gfk = myfprime(x0)
    k = 0
    xk = x0

    # Sets the initial step guess to dx ~ 1
    old_fval = f(xk)
    old_old_fval = old_fval + np.linalg.norm(gfk) / 2

    if retall:
        allvecs = [xk]
    warnflag = 0
    pk = -gfk
    gnorm = vecnorm(gfk, ord=norm)

    sigma_3 = 0.01

    while (gnorm > gtol) and (k < maxiter):
        (dim, mu, shift) = args
        args = (dim, mu * beta, shift)
        # print(k, mu)

        deltak = numpy.dot(gfk, gfk)

        cached_step = [None]

        def polak_ribiere_powell_step(alpha, gfkp1=None):
            xkp1 = xk + alpha * pk
            if gfkp1 is None:
                gfkp1 = myfprime(xkp1)
            yk = gfkp1 - gfk
            beta_k = max(0, numpy.dot(yk, gfkp1) / deltak)
            pkp1 = -gfkp1 + beta_k * pk
            gnorm = vecnorm(gfkp1, ord=norm)
            return (alpha, xkp1, pkp1, gfkp1, gnorm)

        def descent_condition(alpha, xkp1, fp1, gfkp1):
            # Polak-Ribiere+ needs an explicit check of a sufficient
            # descent condition, which is not guaranteed by strong Wolfe.
            #
            # See Gilbert & Nocedal, "Global convergence properties of
            # conjugate gradient methods for optimization",
            # SIAM J. Optimization 2, 21 (1992).
            cached_step[:] = polak_ribiere_powell_step(alpha, gfkp1)
            alpha, xk, pk, gfk, gnorm = cached_step

            # Accept step if it leads to convergence.
            if gnorm <= gtol:
                return True

            # Accept step if sufficient descent condition applies.
            return numpy.dot(pk, gfk) <= -sigma_3 * numpy.dot(gfk, gfk)

        try:
            alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = _line_search_wolfe12(
                f,
                myfprime,
                xk,
                pk,
                gfk,
                old_fval,
                old_old_fval,
                c2=0.4,
                amin=1e-100,
                amax=1e100,
                extra_condition=descent_condition,
            )
        except _LineSearchError:
            # Line search failed to find a better solution.
            warnflag = 2
            break

        # Reuse already computed results if possible
        if alpha_k == cached_step[0]:
            alpha_k, xk, pk, gfk, gnorm = cached_step
        else:
            alpha_k, xk, pk, gfk, gnorm = polak_ribiere_powell_step(alpha_k, gfkp1)

        if retall:
            allvecs.append(xk)
        if callback is not None:
            callback(xk)
        k += 1

    fval = old_fval
    if warnflag == 2:
        msg = _status_message["pr_loss"]
    elif k >= maxiter:
        warnflag = 1
        msg = _status_message["maxiter"]
    else:
        msg = _status_message["success"]

    if disp:
        print("%s%s" % ("Warning: " if warnflag != 0 else "", msg))
        print("         Current function value: %f" % fval)
        print("         Iterations: %d" % k)
        print("         Function evaluations: %d" % func_calls[0])
        print("         Gradient evaluations: %d" % grad_calls[0])

    result = OptimizeResult(
        fun=fval,
        jac=gfk,
        nfev=func_calls[0],
        njev=grad_calls[0],
        status=warnflag,
        success=(warnflag == 0),
        message=msg,
        x=xk,
        nit=k,
    )
    if retall:
        result["allvecs"] = allvecs
    return result

