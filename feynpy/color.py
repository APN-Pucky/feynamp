from re import I

from sympy import Dummy, Function, I, Wild, cacheit, symbols

# from symengine import *


N_c = symbols("N_c", integer=True)
N_g = symbols("N_g", integer=True)
C_A = symbols("C_A", integer=True)
C_F = symbols("C_F", integer=True)
T_F = symbols("T_F", integer=True)

T = Function("T")
f = Function("f")

delta_c = Function("delta_c")
delta_g = Function("delta_g")


def apply_color(expr):
    expr = apply_id(expr)
    expr = apply_f(expr)
    expr = apply_TT(expr)
    expr = apply_delta(expr)
    expr = apply_T(expr)
    expr = apply_id(expr)
    return expr


def apply_id(expr):
    wi = Wild("wi")
    expr = expr.replace(delta_c(wi, wi), N_c)
    expr = expr.replace(delta_g(wi, wi), N_g)
    return expr


def apply_delta(expr, expand=True):
    wa, wb, wc, wd, wi, wj, wk, ww = symbols("wa wb wc wd wi wj wk ww", cls=Wild)
    (
        d1,
        d2,
        d3,
        d4,
        d5,
        d6,
        d7,
        d8,
    ) = symbols("d1 d2 d3 d4 d5 d6 d7 d8", cls=Dummy)
    if expand:
        expr = expr.expand()  # so that we hit it

    expr = expr.replace(delta_c(wa, wb) * delta_c(wb, wc) * ww, delta_c(wa, wc) * ww)
    expr = expr.replace(delta_c(wa, wb) * delta_c(wc, wb) * ww, delta_c(wa, wc) * ww)
    expr = expr.replace(delta_c(wb, wa) * delta_c(wb, wc) * ww, delta_c(wa, wc) * ww)
    expr = expr.replace(delta_c(wb, wa) * delta_c(wc, wb) * ww, delta_c(wa, wc) * ww)

    expr = expr.replace(ww * T(wa, wi, wj) * delta_c(wi, wk), ww * T(wa, wk, wj))
    expr = expr.replace(ww * T(wa, wj, wi) * delta_c(wi, wk), ww * T(wa, wj, wk))
    expr = expr.replace(ww * T(wa, wi, wj) * delta_c(wk, wi), ww * T(wa, wk, wj))
    expr = expr.replace(ww * T(wa, wj, wi) * delta_c(wk, wi), ww * T(wa, wj, wk))

    expr = expr.replace(ww * T(wa, wi, wj) * delta_g(wa, wb), ww * T(wb, wi, wj))
    expr = expr.replace(ww * T(wa, wi, wj) * delta_g(wb, wa), ww * T(wb, wi, wj))

    # expr = expr.replace(ww*T(wa,wi,wj)* delta_c(wi,wj), ww*T(wa,d1,d1))
    # expr = expr.replace(ww*T(wa,wi,wj)* delta_c(wj,wi), ww*T(wa,d2,d2))

    expr = expr.replace(ww * f(wa, wb, wc) * delta_g(wa, wd), ww * f(wd, wb, wc))
    expr = expr.replace(ww * f(wa, wb, wc) * delta_g(wb, wd), ww * f(wa, wd, wc))
    expr = expr.replace(ww * f(wa, wb, wc) * delta_g(wc, wd), ww * f(wa, wb, wd))

    expr = expr.replace(ww * f(wa, wb, wc) * delta_g(wd, wa), ww * f(wd, wb, wc))
    expr = expr.replace(ww * f(wa, wb, wc) * delta_g(wd, wb), ww * f(wa, wd, wc))
    expr = expr.replace(ww * f(wa, wb, wc) * delta_g(wd, wc), ww * f(wa, wb, wd))

    expr = expr.replace(ww * f(wa, wb, wc) * delta_c(wa, wb), ww * f(d3, d3, wc))
    expr = expr.replace(ww * f(wa, wb, wc) * delta_c(wa, wc), ww * f(d4, wb, d4))
    expr = expr.replace(ww * f(wa, wb, wc) * delta_c(wb, wc), ww * f(wa, d5, d5))

    expr = expr.replace(ww * f(wa, wb, wc) * delta_c(wb, wa), ww * f(d6, d6, wc))
    expr = expr.replace(ww * f(wa, wb, wc) * delta_c(wc, wa), ww * f(d7, wb, d7))
    expr = expr.replace(ww * f(wa, wb, wc) * delta_c(wc, wb), ww * f(wa, d8, d8))
    return expr


def apply_T(expr):
    wa, wi = symbols("wa wi", cls=Wild)
    # traceless
    expr = expr.replace(T(wa, wi, wi), 0)
    return expr


def apply_f(expr):
    wa, wb, wc, ww = symbols("wa wb wc ww", cls=Wild)
    i, j, k = symbols("i j k", cls=Dummy)
    dict = expr.match(ww * f(wa, wb, wc))
    expr = expr.replace(
        f(wa, wb, wc) * ww,
        2
        * I
        * (
            T(wc, i, j) * T(wb, j, k) * T(wa, k, i)
            - T(wa, i, j) * T(wb, j, k) * T(wc, k, i)
        )
        * ww,
    )
    return expr


def apply_TT(expr):
    """

    Examples
    ========

    >>> from feynpy.color import *
    >>> i,j,k,l,g = symbols('i j k l g')
    >>> m = T(g,i,j)*T(g,k,l)
    >>> apply_TT(m)
    >>> m
    -C_A*Color.delta(i, k)*Color.delta(j, l)/2 + Color.delta(i, l)*Color.delta(j, k)/2

    """
    wi, wj, wk, wl, wg, ww = symbols("wi wj wk wl wg ww", cls=Wild)
    expr = expr.replace(
        T(wg, wi, wj) * T(wg, wk, wl) * ww,
        (
            delta_c(wi, wl) * delta_c(wk, wj) / 2
            - 1 / N_c / 2 * delta_c(wi, wj) * delta_c(wk, wl)
        )
        * ww,
    )
    return expr
