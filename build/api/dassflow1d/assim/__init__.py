from .m1qn3 import minimize_m1qn3


def minimize(fun, x, jac=False, status=False, tol=1e-6, callback=None, method="M1QN3", **kwargs):
    
    if method == "M1QN3":
        return minimize_m1qn3(fun, x, jac, status, tol, callback, **kwargs)
    else:
        raise ValueError("Unknown 'method' %s" % method)
