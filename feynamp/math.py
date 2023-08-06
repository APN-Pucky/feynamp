from feynamp.leg import get_leg_math
from feynamp.propagator import get_propagator_math
from feynamp.vertex import get_vertex_math


def get_math(object, fd, model):
    if object in fd.vertices:
        return get_vertex_math(object, fd, model)
    elif object in fd.legs:
        return get_leg_math(object, fd, model)
    elif object in fd.propagators:
        return get_propagator_math(object, fd, model)
    else:
        raise Exception("Object not found in feynman diagram")
