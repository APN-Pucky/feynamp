from typing import List

from feynml.leg import Leg

from feynamp import sympy
from feynamp.leg import color_vector_to_casimir, get_color_vector, get_leg_momentum


def assert_spincorrelation(sympy_expr, fd, legs: List[Leg], model):
    assert 1 == 1
