from sympy.parsing.sympy_parser import parse_expr

from feynamp.sympy.color import identity
from feynamp.sympy.lorentz import gamma


def string_to_sympy(s):
    s = s.replace(
        "Gamma", "gamma"
    )  # we keep string like the ufo, but want lowercase gamma for sympy
    s = s.replace("Metric", "metric")
    s = s.replace("complex(0,1)", "I")  # sympy uses I for imaginary unit
    return parse_expr(s, local_dict={"gamma": gamma, "Identity": identity})
    # return parse_expr(s, evaluate=False)


def fancy_print(expr):
    s = string_to_sympy(expr)
    s.replace()
