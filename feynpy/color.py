import itertools

from sympy import Function, Wild, symbols
from sympy.tensor.tensor import (
    TensorHead,
    TensorIndex,
    TensorIndexType,
    WildTensorIndex,
    tensor_heads,
    tensor_indices,
)

N_c = symbols("N_c", integer=True)
N_g = symbols("N_g", integer=True)
C_A = symbols("C_A", integer=True)
C_F = symbols("C_F", integer=True)
T_F = symbols("T_F", integer=True)

ColorIndex = TensorIndexType(
    "Color", dummy_name="C", dim=N_c, metric_symmetry=0
)  # Bug: None crashes replacements
GluonIndex = TensorIndexType(
    "Gluon", dummy_name="G", dim=N_g, metric_symmetry=0
)  # Bug: None crashes replacements

T = TensorHead("T", [GluonIndex, ColorIndex, ColorIndex])
f = TensorHead("f", [GluonIndex, GluonIndex, GluonIndex])

T = Function("T")
f = Function("f")
delta = Function("delta")


def apply_id(expr):
    i = TensorIndex(True, ColorIndex)
    j = TensorIndex(True, GluonIndex)
    expr = expr.subs(ColorIndex.delta(-i, i), ColorIndex.dim)
    expr = expr.subs(ColorIndex.delta(i, -i), ColorIndex.dim)
    expr = expr.subs(GluonIndex.delta(-j, j), GluonIndex.dim)
    expr = expr.subs(GluonIndex.delta(j, -j), GluonIndex.dim)
    return expr


def apply_TT(expr):
    """

    Examples
    ========

    >>> from feynpy.color import *
    >>> i,j,k,l = tensor_indices('i,j,k,l',ColorIndex)
    >>> g = tensor_indices('g',GluonIndex)
    >>> m = T(g,i,j)* T(-g,k,l)
    >>> apply_TT(m)
    >>> m
    -C_A*Color.delta(i, k)*Color.delta(j, l)/2 + Color.delta(i, l)*Color.delta(j, k)/2

    """
    i = WildTensorIndex("i", ColorIndex, ignore_updown=True)
    j = WildTensorIndex("j", ColorIndex, ignore_updown=True)
    k = WildTensorIndex("k", ColorIndex, ignore_updown=True)
    l = WildTensorIndex("l", ColorIndex, ignore_updown=True)
    t = WildTensorIndex("t", GluonIndex, ignore_updown=True)
    rem = Wild("rem")
    expr = expr.replace(
        rem * T(t, i, j) * T(-t, l, k),
        rem
        * (
            ColorIndex.delta(i, l) * ColorIndex.delta(k, j) / 2
            - 1 / N_c / 2 * ColorIndex.delta(i, j) * ColorIndex.delta(k, l)
        ),
    )
    return expr
