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
identity = Function("Identity")

VC = Function("VC")
VA = Function("VA")


def apply_color(expr, expand=True):
    # if expr is tuple
    if isinstance(expr, tuple):
        expr = expr[0]
    expr = apply_color_sum(expr)
    expr = apply_id(expr)
    expr = apply_f(expr)
    expr = apply_TT(expr)
    expr = apply_T(expr)
    expr = apply_id(expr)
    if expand:
        expr = expr.expand()  # so that we hit it
    expr = apply_delta_delta(expr)
    expr = apply_T_delta(expr)
    expr = apply_f_delta(expr)
    expr = apply_T(expr)
    expr = apply_id(expr)
    return expr


def apply_color_sum(expr):
    wi, wj, wk, ww = symbols("wi wj wk ww", cls=Wild)
    expr = expr.replace(VC(wi, wk) * VC(wj, wk) * ww, ww * delta_c(wi, wj))
    expr = expr.replace(VA(wi, wk) * VA(wj, wk) * ww, ww * delta_g(wi, wj))
    return expr


def apply_delta(expr):
    expr = apply_delta_delta(expr)
    expr = apply_T_delta(expr)
    expr = apply_f_delta(expr)

    return expr  # , bool(chagned1 or chagned2 or chagned3)


def apply_id(expr):
    wi = Wild("wi")
    expr, map1 = expr.replace(delta_c(wi, wi), N_c, map=True)
    expr, map2 = expr.replace(delta_g(wi, wi), N_g, map=True)
    return expr  # , bool(map1 or map2)


def apply_delta_delta(expr):
    wa, wb, wc, wd, wi, wj, wk, ww = symbols("wa wb wc wd wi wj wk ww", cls=Wild)
    expr, map1 = expr.replace(
        delta_c(wa, wb) * delta_c(wb, wc) * ww, delta_c(wa, wc) * ww, map=True
    )
    expr, map2 = expr.replace(
        delta_c(wa, wb) * delta_c(wc, wb) * ww, delta_c(wa, wc) * ww, map=True
    )
    expr, map3 = expr.replace(
        delta_c(wb, wa) * delta_c(wb, wc) * ww, delta_c(wa, wc) * ww, map=True
    )
    expr, map4 = expr.replace(
        delta_c(wb, wa) * delta_c(wc, wb) * ww, delta_c(wa, wc) * ww, map=True
    )
    return expr  # , bool(map1 or map2 or map3 or map4)


def apply_T_delta(expr):
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
    expr, map1 = expr.replace(
        ww * T(wa, wj, wi) * delta_c(wi, wk), ww * T(wa, wj, wk), map=True
    )
    expr, map2 = expr.replace(
        ww * T(wa, wi, wj) * delta_c(wi, wk), ww * T(wa, wk, wj), map=True
    )
    expr, map3 = expr.replace(
        ww * T(wa, wi, wj) * delta_c(wk, wi), ww * T(wa, wk, wj), map=True
    )
    expr, map4 = expr.replace(
        ww * T(wa, wj, wi) * delta_c(wk, wi), ww * T(wa, wj, wk), map=True
    )

    expr, map5 = expr.replace(
        ww * T(wa, wi, wj) * delta_g(wa, wb), ww * T(wb, wi, wj), map=True
    )
    expr, map6 = expr.replace(
        ww * T(wa, wi, wj) * delta_g(wb, wa), ww * T(wb, wi, wj), map=True
    )

    # expr = expr.replace(ww*T(wa,wi,wj)* delta_c(wi,wj), ww*T(wa,d1,d1))
    # expr = expr.replace(ww*T(wa,wi,wj)* delta_c(wj,wi), ww*T(wa,d2,d2))
    return expr  # , bool(map1 or map2 or map3 or map4 or map5 or map6)


def apply_f_delta(expr):
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
    expr, map1 = expr.replace(
        ww * f(wa, wb, wc) * delta_g(wa, wd), ww * f(wd, wb, wc), map=True
    )
    expr, map2 = expr.replace(
        ww * f(wa, wb, wc) * delta_g(wb, wd), ww * f(wa, wd, wc), map=True
    )
    expr, map3 = expr.replace(
        ww * f(wa, wb, wc) * delta_g(wc, wd), ww * f(wa, wb, wd), map=True
    )

    expr, map4 = expr.replace(
        ww * f(wa, wb, wc) * delta_g(wd, wa), ww * f(wd, wb, wc), map=True
    )
    expr, map5 = expr.replace(
        ww * f(wa, wb, wc) * delta_g(wd, wb), ww * f(wa, wd, wc), map=True
    )
    expr, map6 = expr.replace(
        ww * f(wa, wb, wc) * delta_g(wd, wc), ww * f(wa, wb, wd), map=True
    )

    expr, map7 = expr.replace(
        ww * f(wa, wb, wc) * delta_c(wa, wb), ww * f(d3, d3, wc), map=True
    )
    expr, map8 = expr.replace(
        ww * f(wa, wb, wc) * delta_c(wa, wc), ww * f(d4, wb, d4), map=True
    )
    expr, map9 = expr.replace(
        ww * f(wa, wb, wc) * delta_c(wb, wc), ww * f(wa, d5, d5), map=True
    )

    expr, map10 = expr.replace(
        ww * f(wa, wb, wc) * delta_c(wb, wa), ww * f(d6, d6, wc), map=True
    )
    expr, map11 = expr.replace(
        ww * f(wa, wb, wc) * delta_c(wc, wa), ww * f(d7, wb, d7), map=True
    )
    expr, map12 = expr.replace(
        ww * f(wa, wb, wc) * delta_c(wc, wb), ww * f(wa, d8, d8), map=True
    )
    return expr  # , bool(map1 or map2 or map3 or map4 or map5 or map6 or map7 or map8 or map9 or map10 or map11 or map12)


def apply_T(expr):
    wa, wi = symbols("wa wi", cls=Wild)
    # traceless
    expr, map = expr.replace(T(wa, wi, wi), 0, map=True)
    return expr  # , bool(map)


def apply_f(expr):
    wa, wb, wc, ww = symbols("wa wb wc ww", cls=Wild)
    i, j, k = symbols("i j k", cls=Dummy)
    expr, map = expr.replace(
        f(wa, wb, wc) * ww,
        2
        * I
        * (
            T(wc, i, j) * T(wb, j, k) * T(wa, k, i)
            - T(wa, i, j) * T(wb, j, k) * T(wc, k, i)
        )
        * ww,
        map=True,
    )
    return expr  # , bool(map)


def apply_TT(expr):
    """

    Examples
    ========

    >>> from feynpy.color import *
    >>> i,j,k,l,g = symbols('i j k l g')
    >>> m = T(g,i,j)*T(g,k,l)
    >>> apply_TT(m)
    >>> m
    -1/C_A*Color.delta(i, k)*Color.delta(j, l)/2 + Color.delta(i, l)*Color.delta(j, k)/2

    """
    wi, wj, wk, wl, wg, ww = symbols("wi wj wk wl wg ww", cls=Wild)
    expr, map = expr.replace(
        T(wg, wi, wj) * T(wg, wk, wl) * ww,
        (
            delta_c(wi, wl) * delta_c(wk, wj) / 2
            - 1 / N_c / 2 * delta_c(wi, wj) * delta_c(wk, wl)
        )
        * ww,
        map=True,
    )
    return expr  # , bool(map)
