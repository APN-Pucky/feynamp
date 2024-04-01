from feynamp.sympy.color import *


def test_apply_TT():
    i, j, k, l, g = symbols("i j k l g")
    m = T(g, i, j) * T(g, k, l)
    result = apply_TT(m)
    expected = delta_c(i, l) * delta_c(k, j) / 2 - delta_c(i, j) * delta_c(k, l) / (
        2 * N_c
    )
    assert expected.equals(result)


def test_apply_TT_scaled():
    i, j, k, l, g = symbols("i j k l g")
    m = T(g, i, j) * T(g, k, l) * 5
    result = apply_TT(m)
    expected = (
        delta_c(i, l) * delta_c(k, j) / 2 * 5
        - delta_c(i, j) * delta_c(k, l) / (2 * N_c) * 5
    )
    assert expected.equals(result)


# TODO temporarily commentent because it is to slow
#
def test_fTT():
    a, b, c, i, j, k, l, g = symbols("a b c i j k l g")
    test = f(a, b, c) * T(c, i, j) * T(b, k, i)
    result = apply_color(apply_color(apply_color(test)))
    expected = I / 2 * N_c * T(a, k, j)
    assert expected.equals(result)


def test_TTT():
    a, b, c, i, j, k, l, g = symbols("a b c i j k l g")
    test = T(b, j, i) * T(b, l, k) * T(a, k, j)
    result = apply_color(apply_color(test))
    expected = -T(a, l, i) / 2 / N_c
    assert expected.equals(result)
