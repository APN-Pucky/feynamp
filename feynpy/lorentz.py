from re import I

from sympy import Dummy, Function, I, Wild, cacheit, symbols

eps = Function("eps")
eps_star = Function("eps_star")
u = Function("u")
u_bar = Function("u_bar")
v = Function("v")
v_bar = Function("v_bar")

gamma = Function("gamma")
gamma_id = Function("gamma_id")
ProjP = Function("P_+")
ProjM = Function("P_-")
metric = Function("eta")

P = Function("P")


def apply_polarisation_sum(expr):
    wmu, wnu, wc, wd, wi, wj, wk, ww = symbols("wmu wnu wc wd wi wj wk ww", cls=Wild)
    expr = expr.replace(
        eps_star(wmu, wd, wc) * eps(wnu, wi, wc), -metric(wmu, wnu), map=True
    )
    # TODO do the pol sums
    return expr


def apply_dirac_trick(expr):
    wa, wb, wc, wd, wi, wj, wk, ww = symbols("wa wb wc wd wi wj wk ww", cls=Wild)
    maps = True
    while maps:
        di, dj = symbols("di dj", cls=Dummy)
        expr, maps = expr.replace(
            u(wb, wa) * u_bar(wc, wa) * ww,
            gamma(di, wb, wc) * P(di, wa) * ww
            + gamma_id(wb, wc) * P(di, wa) * P(di, wa) * ww,
            map=True,
        )
    maps = True
    while maps:
        di, dj = symbols("di dj", cls=Dummy)
        expr, maps = expr.replace(
            v(wb, wa) * v_bar(wc, wa) * ww,
            gamma(dj, wb, wc) * P(dj, wa) * ww
            - gamma_id(wb, wc) * P(dj, wa) * P(dj, wa) * ww,
            map=True,
        )
    return expr


def apply_gammas(expr):
    expr = apply_5_gamma(expr)
    expr = apply_4_gamma(expr)
    expr = apply_3_gamma(expr)
    expr = apply_2_gamma(expr)
    return expr


def apply_2_gamma(expr):
    wa, wb, wc, wd, wi, wj, wk, ww = symbols("wa wb wc wd wi wj wk ww", cls=Wild)
    maps = True
    while maps:
        di, dj = symbols("di dj", cls=Dummy)
        expr, maps = expr.replace(
            gamma(wa, wb, wd) * gamma(wa, wd, wc) * ww, -gamma_id(wb, wc) * ww, map=True
        )
    return expr


def apply_3_gamma(expr):
    wa, wb, wc, wd, we, wg, wi, wj, wk, ww = symbols(
        "wa wb wc wd we wg wi wj wk ww", cls=Wild
    )
    maps = True
    while maps:
        di, dj = symbols("di dj", cls=Dummy)
        expr, maps = expr.replace(
            gamma(wa, wb, wd) * gamma(wg, wd, wc) * gamma(wa, wc, we) * ww,
            -2 * gamma(wg, wb, we) * ww,
            map=True,
        )
    return expr


def apply_4_gamma(expr):
    wa, wb, wc, wd, we, wf, wg, wh, wi, wj, wk, ww = symbols(
        "wa wb wc wd we wf wg wh wi wj wk ww", cls=Wild
    )
    maps = True
    while maps:
        di, dj = symbols("di dj", cls=Dummy)
        expr, maps = expr.replace(
            gamma(wa, wb, wd)
            * gamma(wi, wd, we)
            * gamma(wj, we, wf)
            * gamma(wa, wf, wc)
            * ww,
            4 * metric(wi, wj) * gamma_id(wb, wc) * ww,
            map=True,
        )
    return expr


def apply_5_gamma(expr):
    wa, wb, wc, wd, we, wf, wg, wh, wi, wj, wk, ww = symbols(
        "wa wb wc wd we wf wg wh wi wj wk ww", cls=Wild
    )
    maps = True
    while maps:
        print(maps)
        di, dj = symbols("di dj", cls=Dummy)
        expr, maps = expr.replace(
            gamma(wa, wb, wd)
            * gamma(wi, wd, we)
            * gamma(wj, we, wf)
            * gamma(wk, wf, wg)
            * gamma(wa, wg, wh)
            * ww,
            -2 * gamma(wi, wb, di) * gamma(wj, di, dj) * gamma(wk, dj, wh) * ww,
            map=True,
        )
    return expr
