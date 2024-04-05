from .amplitude import complex_conjugate, feynman_diagram_to_string
from .log import debug
from .propagator import find_propagator_in_model, get_propagator_math
from .vertex import find_vertex_in_model, get_vertex_math


def get_spin_average(fds):
    """
    Initial fermions produce an averaging by 1/2 each
    """
    r = []
    for leg in fds[0].get_incoming():
        if leg.is_any_fermion():
            r += ["1/2"]
        elif leg.pdgid == 21:
            r += ["1/2"]
        elif leg.pdgid == 22:
            r += ["1/2"]
        else:
            raise ValueError(f"TODO: Unknown color average for pdgid {leg.pdgid}")
    debug(f"get_spin_average():{r}")
    return r


def get_color_average(fds):
    """
    Initial quarks produce an averaging by 1/3 each
    Initial gluons produce an averaging by 1/8 each
    """
    r = []
    for leg in fds[0].get_incoming():
        if leg.pdgid == 21:
            r += ["1/8"]
        elif leg.pdgid in range(1, 7) or -leg.pdgid in range(1, 7):
            r += ["1/3"]
        elif leg.pdgid == 22:
            r += ["1"]
        elif leg.pdgid == 23:
            r += ["1"]
        elif leg.pdgid in range(11, 19) or -leg.pdgid in range(11, 19):
            r += ["1"]
        else:
            raise ValueError(f"TODO: Unknown color average for pdgid {leg.pdgid}")
    debug(f"get_color_average():{r}")
    return r
