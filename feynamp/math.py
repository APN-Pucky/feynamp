from feynamp.leg import get_leg_math
from feynamp.propagator import get_propagator_math
from feynamp.vertex import get_vertex_math


def get_math(obj, fd, model):
    if obj in fd.vertices:
        return get_vertex_math(obj, fd, model)
    elif obj in fd.legs:
        return get_leg_math(obj, fd, model)
    elif obj in fd.propagators:
        return get_propagator_math(obj, fd, model)
    else:
        raise Exception(f"Object {obj} not found in feynman diagram")
