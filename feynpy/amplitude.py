from feynpy.leg import get_leg_math_string
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
