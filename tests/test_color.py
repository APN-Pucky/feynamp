from feynpy.color import *


def test_apply_TT():
    i, j, k, l = tensor_indices("i,j,k,l", ColorIndex)
    g = tensor_indices("g", GluonIndex)
    m = 5 * T(g, i, j) * T(-g, k, l)
    l = apply_TT(m)
    assert l == 5
