import re

from feynamp.log import warning
from feynamp.momentum import insert_momentum
from feynamp.util import find_particle_in_model


def get_leg_math_string(leg, fd, model):
    return get_leg_math(leg, fd, model)


def get_leg_momentum(leg):
    if leg is None:
        return None
    if leg.momentum is None or leg.momentum.name is None:
        raise ValueError("Momentum not set for particle")
    mom = insert_momentum(leg.momentum.name)
    return mom


def color_vector_to_operator(color_vector):
    if color_vector == "VA":
        return "f"
    if color_vector == "VC":
        return "T"
    return None


def color_vector_to_id(color_vector):
    if color_vector == "VA":
        return "da"
    if color_vector == "VC":
        return "df"
    return None


def color_vector_to_index(color_vector):
    if color_vector == "VA":
        return "Glu"
    if color_vector == "VC":
        return "Color"
    return None


def is_swapped_color_vector(leg, s2):
    """
    For colorcorrelations the chains of T's must be properly ordered so we check if they must be swapped

    TODO: check this for gluons, maybe sign is wrong...
    """
    if re.search(r"T\(Color" + leg.id + r",.*?,.*?\)", s2):
        return True
    elif re.search(r"T\(.*?,Color" + leg.id + r",.*?\)", s2):
        return False
    elif re.search(r"f\(Glu" + leg.id + r",.*?,.*?\)", s2):
        return True
    elif re.search(r"f\(.*?,Glu" + leg.id + r",.*?\)", s2):
        return False
    elif re.search(r"f\(.*?,.*?,Glu" + leg.id + r"\)", s2):
        warning("leg color third in f, check colorcorrelations")
        return True
    raise ValueError(f"Color vector for {leg} not found in squared amplitude")


def get_color_vector(fd, leg, model):
    p = find_leg_in_model(fd, leg, model)
    # give particles color vectors to sum over them in the end (or better average)
    # TODO this could be also done as incoming vs outcoming
    if p.color == 8:
        # if particle is a gluon give it a adjoint color function
        return "VA"
    if p.color == 3 or p.color == -3:
        # if particle is a quark give it a fundamental color function
        return "VC"
    return None


def get_leg_math(fd, leg, model):  # epsilons or u/v optionally also barred
    p = find_leg_in_model(fd, leg, model)
    mom = get_leg_momentum(leg)
    ret = ""
    # give particles color vectors to sum over them in the end (or better average)
    # TODO this could be also done as incoming vs outcoming
    if p.color == 8:
        # if particle is a gluon give it a adjoint color function
        ret += f"VA(Glu{p.particle.id},{mom})*"
    if p.color == 3 or p.color == -3:
        # if particle is a quark give it a fundamental color function
        ret += f"VC(Color{p.particle.id},{mom})*"

    if p.spin == 3:
        if leg.is_incoming():
            ret += f"eps(Mu{p.particle.id},Pol{p.particle.id},{mom})"
        else:
            ret += f"eps_star(Mu{p.particle.id},Pol{p.particle.id},{mom})"
    if p.spin == 2:
        if not p.particle.is_anti():
            if leg.is_incoming():
                ret += f"u(Spin{p.particle.id},{mom})"
            else:
                ret += f"u_bar(Spin{p.particle.id},{mom})"
        else:
            if leg.is_incoming():
                ret += f"v(Spin{p.particle.id},{mom})"
            else:
                ret += f"v_bar(Spin{p.particle.id},{mom})"
    return ret


def find_leg_in_model(fd, leg, model):  # find leg in model
    assert leg in fd.legs
    return find_particle_in_model(leg, model)
