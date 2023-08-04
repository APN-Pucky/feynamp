from sympy.parsing.sympy_parser import parse_expr

from feynpy.leg import get_leg_math_string
from feynpy.lorentz import gamma
from feynpy.propagator import get_propagator_math_string
from feynpy.vertex import get_vertex_math_string


def feynman_diagram_to_string(feynman_diagram, feyn_model):
    fd = feynman_diagram
    vm = []
    lm = []
    pm = []
    for v in fd.vertices:
        vm.append(get_vertex_math_string(fd, v, feyn_model))
    for l in fd.legs:
        lm.append(get_leg_math_string(fd, l, feyn_model))
    for p in fd.propagators:
        pm.append(get_propagator_math_string(fd, p, feyn_model))
    return f"{' * '.join(vm)} * {' * '.join(lm)} * {' * '.join(pm)}"


def string_to_sympy(s):
    s = s.replace(
        "Gamma", "gamma"
    )  # we keep string like the ufo, but want lowercase gamma for sympy
    s = s.replace("Metric", "metric")
    s = s.replace("complex(0,1)", "I")  # sympy uses I for imaginary unit
    return parse_expr(s, local_dict={"gamma": gamma})
    # return parse_expr(s, evaluate=False)
